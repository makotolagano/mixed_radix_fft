"""5G NR EVM test for the fixed-point DFT model (transform-precoder role).

The supported sizes (12 * 2^a * 3^b * 5^c) are exactly the NR DFT-s-OFDM
allocation sizes (TS 38.211 transform precoding), so this sweeps EVERY legal
allocation x modulation and checks the DFT's error contribution against an
allocated slice of the TS 38.101-1 transmitter EVM budget (Table 6.4.2.1-1:
QPSK 17.5%, 16QAM 12.5%, 64QAM 8%, 256QAM 3.5% -- the WHOLE transmitter;
EVM contributions add in power, so a datapath block must stay far below).

Method (component-level 38.101 Annex-B style, ideal channel):
  * stimulus: unit-average-power constellation symbols, scaled so the PEAK
    amplitude sits at --peak (default 0.9) of full scale (input contract)
  * DUT: the model chain exactly as chain_ref.py composes it (ascending
    radices, octant twiddles, the full RTL convention) -- K frames streamed
    back-to-back plus a zero flush frame
  * reference: ideal numpy DFT of the transmitted (float) symbols; the chain
    output is digit-reverse reordered, then a SINGLE complex gain alpha is
    least-squares fitted per frame (3GPP amplitude/phase normalization; a
    flat alpha is stricter than the spec's per-subcarrier equalization and
    absorbs the chain's 2^-total_shift scaling)
  * EVM_rms = sqrt(sum|meas/alpha - ref|^2 / sum|ref|^2) per frame

Run:  .venv/bin/python model/evm_5g.py                (full 53 x 4 sweep)
      .venv/bin/python model/evm_5g.py --sizes 12 3240 --mods qam256
"""
import argparse
import csv
import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from mixed_radix_fft_fxp import MixedRadix_SDF_stage_counter_ctrl_FXP, supported_fft_sizes
from chain_ref import factorize_ascending
from utils import digit_reverse

# TS 38.101-1 Table 6.4.2.1-1 (total transmitter EVM, %)
NR_EVM_LIMITS = {'qpsk': 17.5, 'qam16': 12.5, 'qam64': 8.0, 'qam256': 3.5}

DATA_DTYPE    = 'fxp-s18/16'
COEFF_DTYPE   = 'fxp-s18/16'
TWIDDLE_DTYPE = 'fxp-s18/16'


def constellation(mod):
    """Unit-average-power constellation points (one axis) and peak |symbol|."""
    levels = {'qpsk': 1, 'qam16': 2, 'qam64': 3, 'qam256': 4}[mod]
    pts = np.arange(-(2 ** levels - 1), 2 ** levels, 2, dtype=float)
    pts = pts / np.sqrt(np.mean(pts ** 2) * 2)      # unit average power (I+Q)
    peak = np.sqrt(2.0) * np.max(np.abs(pts))
    return pts, peak


def gen_symbols(mod, n, rng):
    pts, _ = constellation(mod)
    return pts[rng.integers(0, len(pts), n)] + 1j * pts[rng.integers(0, len(pts), n)]


def build_chain(n):
    radices = factorize_ascending(n)
    sizes, acc = [], 1
    for r in radices:
        sizes.append(n // acc)
        acc *= r
    stages = [
        MixedRadix_SDF_stage_counter_ctrl_FXP(
            config=r, stage_index=0, size=s, cfg_delay=s // r,
            dtype=DATA_DTYPE, coeff_dtype=COEFF_DTYPE, preadder_shift_bits=None,
            overflow='wrap', rounding='half_up', capability=r,
            twiddle_source='octant', twiddle_dtype=TWIDDLE_DTYPE,
            preadder_exit_round=True, preadder_internal_frac=22)
        for r, s in zip(radices, sizes)
    ]
    lat = sum((r - 1) * (s // r) for r, s in zip(radices, sizes))
    rev = digit_reverse(list(reversed(radices)))
    return stages, lat, rev


def evm_point(n, mod, frames, peak_level, seed):
    """EVM (%) per frame for one (size, modulation) point."""
    rng = np.random.default_rng(seed + n)
    _, peak = constellation(mod)
    scale = peak_level / peak

    x = np.concatenate([gen_symbols(mod, n, rng) for _ in range(frames)]) * scale
    stream = np.concatenate([x, np.zeros(n, dtype=complex)])   # + flush frame

    stages, lat, rev = build_chain(n)
    out = np.empty(len(stream), dtype=complex)
    for i, s_in in enumerate(stream):
        y = s_in
        for st in stages:
            y = st.calculate(y)
        out[i] = y
    out = out[lat:lat + frames * n]

    evms = []
    for f in range(frames):
        ref  = np.fft.fft(x[f * n:(f + 1) * n])
        meas = out[f * n:(f + 1) * n][rev]
        alpha = np.vdot(ref, meas) / np.vdot(ref, ref)   # LS complex gain fit
        err = meas / alpha - ref
        evms.append(100.0 * np.sqrt(np.sum(np.abs(err) ** 2) / np.sum(np.abs(ref) ** 2)))
    return evms


def sweep_size(n, mods, frames, peak_level, seed):
    rows = []
    for mod in mods:
        evms = evm_point(n, mod, frames, peak_level, seed)
        rows.append({'size': n, 'mod': mod, 'frames': frames,
                     'evm_rms_pct': float(np.sqrt(np.mean(np.square(evms)))),
                     'evm_max_pct': float(np.max(evms)),
                     'nr_limit_pct': NR_EVM_LIMITS[mod]})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sizes', type=int, nargs='+', default=None,
                    help='allocation sizes (default: all supported = all legal '
                         'NR transform-precoding sizes)')
    ap.add_argument('--mods', nargs='+', default=list(NR_EVM_LIMITS),
                    choices=list(NR_EVM_LIMITS))
    ap.add_argument('--frames', type=int, default=4, help='frames (seeds) per point')
    ap.add_argument('--peak', type=float, default=0.9,
                    help='constellation peak as fraction of full scale (input contract)')
    ap.add_argument('--budget', type=float, default=1.0,
                    help='allocated EVM budget for this block, %% (PASS criterion)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--jobs', type=int, default=0,
                    help='parallel workers; 0 = one per size capped at half the CPUs')
    ap.add_argument('--csv', default=str(Path(__file__).parent / 'results' / 'evm_5g.csv'))
    a = ap.parse_args()

    sizes = a.sizes or supported_fft_sizes()
    worker = partial(sweep_size, mods=a.mods, frames=a.frames,
                     peak_level=a.peak, seed=a.seed)
    jobs = a.jobs or min(len(sizes), max(1, (os.cpu_count() or 2) // 2))

    rows = []
    if jobs > 1:
        with Pool(processes=jobs) as pool:
            for r in pool.imap(worker, sizes):
                rows.extend(r)
                print(f"N={r[0]['size']:5d}: " + '  '.join(
                    f"{x['mod']}={x['evm_rms_pct']:.4f}%" for x in r), flush=True)
    else:
        for n in sizes:
            r = worker(n)
            rows.extend(r)
            print(f"N={r[0]['size']:5d}: " + '  '.join(
                f"{x['mod']}={x['evm_rms_pct']:.4f}%" for x in r), flush=True)

    os.makedirs(os.path.dirname(a.csv), exist_ok=True)
    with open(a.csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {a.csv}")
    print(f"{'mod':>7} | {'worst EVM (size)':>20} | {'rms EVM':>9} | "
          f"{'NR limit':>8} | {'margin':>8} | vs {a.budget}% budget")
    ok = True
    for mod in a.mods:
        mrows = [r for r in rows if r['mod'] == mod]
        worst = max(mrows, key=lambda r: r['evm_max_pct'])
        rms = float(np.sqrt(np.mean([r['evm_rms_pct'] ** 2 for r in mrows])))
        limit = NR_EVM_LIMITS[mod]
        margin_db = 20 * np.log10(limit / worst['evm_max_pct'])
        verdict = 'PASS' if worst['evm_max_pct'] <= a.budget else 'FAIL'
        ok = ok and (verdict == 'PASS')
        print(f"{mod:>7} | {worst['evm_max_pct']:8.4f}% (N={worst['size']:4d}) | "
              f"{rms:8.4f}% | {limit:7.1f}% | {margin_db:6.1f} dB | {verdict}")
    print("\nEVM " + ("PASS" if ok else "FAIL") +
          f" -- worst-case DFT EVM vs the {a.budget}% block allocation "
          f"(NR limits are for the WHOLE transmitter)")
    raise SystemExit(0 if ok else 1)


if __name__ == '__main__':
    main()
