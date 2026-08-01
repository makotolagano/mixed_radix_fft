"""PYNQ bring-up helpers for the mr_fft IP behind an AXI DMA (64-bit stream).

SELF-CONTAINED: copy this single file to the board (needs numpy + matplotlib
only). The optional bit-exact model check activates automatically when the
repo's model is importable (i.e. on the dev machine).

Stream packing (default): {'0, im, re} -> re in bits [17:0], im in bits
[35:18], zeros above -- both 18-bit two's complement s2.16 codes. If your RTL
slices the 64-bit word differently (e.g. 32-bit aligned fields), change
IM_SHIFT below and nothing else.

The IP output needs one post-step before it looks like an FFT:
  * REORDER: the chain delivers digit-reversed bin order (mixed-radix DIF);
    natural_order() applies the permutation for the configured N
  * (scaling is on-chip: the final scaler delivers the classical X/N
    convention, so reference() is simply fft(x)/N)

Typical notebook flow (configuration + DMA are yours):

    import numpy as np, pynq_fft_test as ft
    from pynq import Overlay, allocate

    ft.list_configs()                            # CONFIG_SEL <-> N table
    n = 300
    sel = ft.config_sel(n)                       # value for the CONFIG_SEL reg
    x = ft.make_signal('tone', n)                # complex float, |peak|<=0.9
    inbuf  = allocate(shape=(n,), dtype=np.uint64)
    outbuf = allocate(shape=(n,), dtype=np.uint64)
    inbuf[:] = ft.pack(x)

    # ... write sel to CONFIG_SEL, CTRL.COMMIT, poll STATUS.BUSY==0 ...
    # ... dma.sendchannel.transfer(inbuf); dma.recvchannel.transfer(outbuf) ...

    meas, ref = ft.analyze(x, outbuf, n)         # plots + SQNR summary

Remember the IP contract: whole frames only, and reconfigure (COMMIT) only
after the previous config has fully drained.
"""
import math

import numpy as np

try:
    import matplotlib.pyplot as plt              # headless-safe: mpl falls back to Agg
except Exception:                                 # plotting optional
    plt = None

# ---------------------------------------------------------------------------
# fixed-point word / stream layout
# ---------------------------------------------------------------------------
WORD_BITS = 18          # s2.16 data word
FRAC_BITS = 16
IM_SHIFT  = 18          # bit position of im within the 64-bit beat ({'0,im,re})

_CODE_MIN = -(1 << (WORD_BITS - 1))
_CODE_MAX = (1 << (WORD_BITS - 1)) - 1
_MASK     = (1 << WORD_BITS) - 1


def quantize(x):
    """Complex float -> integer s2.16 codes (round-to-nearest, saturating)."""
    re = np.clip(np.round(np.real(x) * 2.0 ** FRAC_BITS), _CODE_MIN, _CODE_MAX)
    im = np.clip(np.round(np.imag(x) * 2.0 ** FRAC_BITS), _CODE_MIN, _CODE_MAX)
    return re.astype(np.int64), im.astype(np.int64)


def pack(x):
    """Complex float samples -> uint64 stream beats ({'0, im, re})."""
    re, im = quantize(x)
    return (((im & _MASK) << IM_SHIFT) | (re & _MASK)).astype(np.uint64)


def unpack(beats):
    """uint64 stream beats -> complex float samples (sign-extended codes)."""
    beats = np.asarray(beats, dtype=np.uint64)
    re = (beats & _MASK).astype(np.int64)
    im = ((beats >> IM_SHIFT) & _MASK).astype(np.int64)
    sign = 1 << (WORD_BITS - 1)
    re = (re ^ sign) - sign
    im = (im ^ sign) - sign
    return (re + 1j * im) * 2.0 ** -FRAC_BITS


# ---------------------------------------------------------------------------
# chain output conventions: digit-reversed order + fixed scaling schedule
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# configuration mapping: FFT length N <-> CONFIG_SEL register value
# ---------------------------------------------------------------------------
FFT_SIZE_BASE = 12          # c_fft_size_base in mr_fft_pkg.vhd
MAX_FFT_SIZE  = 3300        # c_max_fft_size


def supported_sizes():
    """All supported FFT lengths, ascending -- IDENTICAL to the RTL's
    c_fft_sizes table (f_fft_sizes sorts ascending: config index = rank of N),
    so list.index(N) IS the CONFIG_SEL register value."""
    sizes, n2 = [], FFT_SIZE_BASE
    while n2 <= MAX_FFT_SIZE:
        n3 = n2
        while n3 <= MAX_FFT_SIZE:
            n5 = n3
            while n5 <= MAX_FFT_SIZE:
                sizes.append(n5)
                n5 *= 5
            n3 *= 3
        n2 *= 2
    return sorted(sizes)


def config_sel(n):
    """CONFIG_SEL register value for FFT length `n` (raises with the list of
    valid lengths if `n` is not a supported configuration)."""
    sizes = supported_sizes()
    try:
        return sizes.index(n)
    except ValueError:
        raise ValueError(f'N={n} is not a supported FFT length; valid: {sizes}') from None


def list_configs():
    """Print the CONFIG_SEL <-> N table with each config's radix factorization
    and scaling (output = X * 2**-shift)."""
    sizes = supported_sizes()
    print(f'{len(sizes)} configurations (CONFIG_SEL -> N):')
    print(f'{"sel":>4} | {"N":>5} | {"radices":<22} | {"stages":>6} | {"2^-shift":>8}')
    print('-' * 60)
    for sel, n in enumerate(sizes):
        radices = factorize(n)
        print(f'{sel:>4} | {n:>5} | {"*".join(map(str, radices)):<22} | '
              f'{len(radices):>6} | 2^-{total_shift(n):<5}')
    return sizes


def factorize(n):
    """Ascending radix list [2..][3..][5..] -- the chain's processing order."""
    radices, rem = [], n
    for r in (2, 3, 5):
        while rem % r == 0:
            radices.append(r)
            rem //= r
    if rem != 1:
        raise ValueError(f'N={n} is not 2^a*3^b*5^c')
    return radices


def total_shift(n):
    """Sum of per-stage scaling shifts (fixed schedule: 2->1, 3->2, 5->3)."""
    return sum({2: 1, 3: 2, 5: 3}[r] for r in factorize(n))


def digit_reverse_perm(radices):
    """Mixed-radix digit-reverse permutation (clone of model/utils.py)."""
    radices = list(radices)
    n = int(np.prod(radices))
    mult, tmp = [], n
    for r in radices:
        tmp //= r
        mult.append(tmp)
    mult_rev, tmp = [], n
    for r in radices[::-1]:
        tmp //= r
        mult_rev.append(tmp)
    perm = np.empty(n, dtype=np.int64)
    for idx in range(n):
        digits, t = [], idx
        for m in mult:
            digits.append(t // m)
            t %= m
        perm[idx] = sum(m * d for m, d in zip(mult_rev, digits[::-1]))
    return perm


def natural_order(out_frame, n):
    """Reorder one captured output frame into natural bin order."""
    rev = digit_reverse_perm(list(reversed(factorize(n))))
    return np.asarray(out_frame)[rev]


def reference(x, n):
    """Ideal spectrum at the IP's scale. The on-chip final scaler
    (mr_fft_scaler) delivers the classical X/N convention."""
    return np.fft.fft(np.asarray(x)) / n


# ---------------------------------------------------------------------------
# test signals (peak <= 0.9 full scale, the chain input contract)
# ---------------------------------------------------------------------------
def make_signal(kind, n, seed=1):
    rng = np.random.default_rng(seed)
    if kind == 'tone':
        # coherent (bin-centered) complex tone at a bin coprime with N: zero
        # leakage, generic twiddle walk -- the spectrum is one clean spike
        k = max(1, n // 7)
        while math.gcd(k, n) != 1:
            k += 1
        idx = np.arange(n)
        return np.sqrt(2) * np.exp(2j * np.pi * k * idx / n)
    if kind == 'qam16':
        pts = np.array([-3, -1, 1, 3], dtype=float) / np.sqrt(10)
        x = pts[rng.integers(0, 4, n)] + 1j * pts[rng.integers(0, 4, n)]
        return x * (0.9 / (np.sqrt(2) * 3 / np.sqrt(10)))   # peak at 0.9
    if kind == 'ramp':
        idx = np.arange(n)
        return ((idx + 1) * 64 - 1j * (idx + 1) * 32) * 2.0 ** -FRAC_BITS
    if kind == 'impulse':
        x = np.zeros(n, dtype=complex)
        x[0] = 0.9
        return x                                            # flat spectrum
    if kind == 'noise':
        return 0.3 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    raise ValueError(f'unknown signal kind: {kind}')


# ---------------------------------------------------------------------------
# analysis / plotting
# ---------------------------------------------------------------------------
def analyze(x, captured, n, save=None, title=None):
    """Unpack + reorder one captured frame, plot it against the ideal
    reference and print an SQNR summary. `captured` may be uint64 beats or
    already-unpacked complex samples. Returns (measured_natural, reference)."""
    captured = np.asarray(captured)
    if captured.dtype == np.uint64:
        captured = unpack(captured[:n])
    meas = natural_order(captured[:n], n)
    ref = reference(x, n)

    err = meas - ref
    p_ref = np.sum(np.abs(ref) ** 2)
    p_err = np.sum(np.abs(err) ** 2)
    sqnr = 10 * np.log10(p_ref / p_err) if p_err > 0 else np.inf
    print(f'N={n}: SQNR vs ideal = {sqnr:.2f} dB '
          f'(EVM ~ {100 * 10 ** (-sqnr / 20):.4f}%)')

    if plt is not None:
        eps = 1e-9
        bins = np.arange(n)
        fig, ax = plt.subplots(3, 1, figsize=(12, 9))

        ax[0].plot(bins, np.real(x), color='#0072B2', linewidth=1.2, label='re')
        ax[0].plot(bins, np.imag(x), color='#E69F00', linewidth=1.2, label='im')
        ax[0].set_title(title or f'mr_fft N={n}')
        ax[0].set_ylabel('input (full scale)')
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(loc='best')

        ax[1].plot(bins, 20 * np.log10(np.abs(ref) + eps), color='black',
                   linewidth=1.6, label='ideal (scaled)')
        ax[1].plot(bins, 20 * np.log10(np.abs(meas) + eps), color='#0072B2',
                   linewidth=1.0, label='IP output (reordered)')
        ax[1].set_ylabel('|X| (dB)')
        ax[1].grid(True, alpha=0.3)
        ax[1].legend(loc='best')

        ax[2].plot(bins, 20 * np.log10(np.abs(err) + eps), color='#D55E00',
                   linewidth=1.0, label='|error|')
        ax[2].set_xlabel('bin (natural order)')
        ax[2].set_ylabel('error (dB)')
        ax[2].ticklabel_format(useOffset=False, axis='y')
        ax[2].grid(True, alpha=0.3)
        ax[2].legend(loc='best')

        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=160)
            print(f'plot -> {save}')
        else:
            plt.show()
        plt.close(fig)

    return meas, ref


# ---------------------------------------------------------------------------
# self-test / dev-machine demo (no hardware): pack/unpack roundtrip, and --
# when the repo model is importable -- a true fixed-point chain run
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # pack/unpack roundtrip over random codes
    rng = np.random.default_rng(1)
    codes = rng.integers(_CODE_MIN, _CODE_MAX + 1, 1000)
    x = (codes + 1j * rng.integers(_CODE_MIN, _CODE_MAX + 1, 1000)) * 2.0 ** -FRAC_BITS
    assert np.array_equal(unpack(pack(x)), x), 'pack/unpack roundtrip failed'
    print('pack/unpack roundtrip: OK (1000 random full-range codes)')

    n = 300
    sig = make_signal('tone', n)
    try:
        from chain_ref import run_cfg
        w = WORD_BITS
        cin = [(int(r), int(i)) for r, i in zip(*quantize(sig))]
        exp = run_cfg(n, cin + [(0, 0)] * n, f'fxp-s{w}/{FRAC_BITS}',
                      'fxp-s18/16', 'fxp-s18/16', FRAC_BITS)[:n]
        captured = np.array([re + 1j * im for re, im in exp]) * 2.0 ** -FRAC_BITS
        print('model chain available: analyzing a REAL fixed-point run')
    except ImportError:
        # simulate a captured (digit-reversed) frame from the ideal spectrum
        perm = digit_reverse_perm(list(reversed(factorize(n))))
        captured = np.empty(n, dtype=complex)
        captured[perm] = reference(sig, n)
        print('model not importable: demo with ideal data')

    analyze(sig, captured, n, save='pynq_fft_test_demo.png')
