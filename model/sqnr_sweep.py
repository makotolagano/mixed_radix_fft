"""SQNR regression sweep: FFT size x datapath rounding mode.

For every requested size N the script factorizes N into an ascending 2/3/5
stage chain (the hardware order), runs the fixed-point chain model at the RTL conventions
(data s18/16, coeff/twiddle s18/16, exit rounding with internal_frac 22,
fixed scaling schedule -- docs/datapath_width_convention.md) once per
rounding mode, and reports SQNR against the NumPy reference. Twiddle ROM
values always stay round-to-nearest (design-time constants), matching the
RTL -- only the datapath/coefficient rounding is swept.

Outputs (under model/results/):
  sqnr_rounding_sweep.csv   one row per (size, mode)
  sqnr_rounding_sweep.png   SQNR vs N, one line per mode

Run:  .venv/bin/python model/sqnr_sweep.py                     (quick set)
      .venv/bin/python model/sqnr_sweep.py --all               (all 53 supported sizes, slow)
      .venv/bin/python model/sqnr_sweep.py --sizes 60 360 3240 --signal qam16
"""
import argparse
import csv
import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.fft import fft

from mixed_radix_fft_fxp import supported_fft_sizes
from test_sqnr_plot import (run_chain, generate_signal, shift_schedule,
                            total_chain_latency, _stage_sizes_from_radices, plt)
from utils import digit_reverse

# spans small->large and pure-2/3/5-heavy mixes; all are supported 12*p sizes
QUICK_SIZES = [12, 24, 48, 60, 96, 144, 240, 360, 540, 900, 1200, 1440, 1536, 2160, 3240]

ROUNDING_MODES = ['floor', 'around', 'half_up']
# fixed categorical order (Okabe-Ito, CVD-safe); markers as secondary encoding
MODE_STYLE = {
    'floor':   ('#0072B2', 'o'),
    'around':  ('#E69F00', 's'),
    'half_up': ('#009E73', '^'),
}


def factorize_235(n):
    """n as an ASCENDING radix list [2..,3..,5..] -- the hardware chain's
    processing order (mid-bypass slot layout; same as chain_ref/pynq/evm
    flows). Stage order matters: descending measures ~13 dB better on tones
    at N=3240, so a descending sweep is not hardware-representative."""
    radices, rem = [], n
    for r in (2, 3, 5):
        while rem % r == 0:
            radices.append(r)
            rem //= r
    if rem != 1:
        raise ValueError(f'N={n} is not 2^a*3^b*5^c')
    return radices


# signals that depend on np.random -- these get --frames seed-averaging;
# deterministic signals always run one frame
RANDOM_SIGNALS = {'qam16'}


def sweep_size(n, signal_kind, modes, dtype, twiddle_dtype, coeff_dtype, seed,
               scaling='fixed', preadder_round='node', internal_frac=None, frames=1):
    radices = factorize_235(n)
    sizes = _stage_sizes_from_radices(radices)
    lat = total_chain_latency(radices, sizes)
    rev = digit_reverse(list(reversed(radices)))
    lsb = 2.0 ** -int(dtype.split('/')[1])
    n_frames = frames if signal_kind in RANDOM_SIGNALS else 1
    # HARDWARE convention: no output rescaler exists in the RTL, so measure
    # the raw chain output against fft(x) * 2**(-total_shift) -- this makes
    # the sweep reproduce PYNQ board measurements exactly
    ref_scale = 2.0 ** -sum(shift_schedule(radices, scaling))

    rows = []
    for mode in modes:
        # power-domain accumulation across frames (proper multi-seed average)
        p_sig = p_noise = 0.0
        bias = 0.0 + 0.0j
        peaks = []
        tone_p_sig = tone_p_spur = 0.0
        for frame in range(n_frames):
            np.random.seed(seed + frame)
            x = generate_signal(signal_kind, n)
            np_ref = fft(x) * ref_scale
            _, fxp_out, widths = run_chain(config=radices, stage_sizes=sizes,
                                           input_signal=np.append(x, np.zeros(n)),
                                           dtype=dtype, twiddle_dtype=twiddle_dtype,
                                           coeff_dtype=coeff_dtype, rounding=mode,
                                           scaling=scaling,
                                           preadder_exit_round=(preadder_round == 'exit'),
                                           preadder_internal_frac=internal_frac,
                                           run_fp=False, final_scaler=False)
            fxp_ss = fxp_out[lat:lat + n][rev]
            err = np_ref - fxp_ss
            p_sig += float(np.sum(np.abs(np_ref) ** 2))
            p_noise += float(np.sum(np.abs(err) ** 2))
            bias += np.mean(err)
            peaks += [v for key, v in widths.items() if key.endswith('_peak')]
            if signal_kind == 'single_tone':
                # For a coherent DIGITAL tone the reference has zero leakage, so
                # sqnr_db already IS the datasheet tone SNR (tone power vs total
                # error power -- including the fundamental bin, where this SDF
                # chain deposits nearly all its rounding error as gain/phase
                # error; the off-bin residues largely cancel exactly). SFDR
                # (tone vs largest single off-bin spur) is kept as the spur
                # diagnostic; it can be inf when every off-bin rounds to zero.
                kbin = int(np.argmax(np.abs(np_ref)))
                tone_p_sig += float(np.abs(np_ref[kbin]) ** 2)
                tone_p_spur = max(tone_p_spur,
                                  float(np.max(np.abs(np.delete(fxp_ss, kbin)) ** 2)))

        def db(num, den):
            return float('inf') if den == 0 else 10 * np.log10(num / den)

        bias = bias / n_frames / lsb
        rows.append({'size': n, 'radices': '*'.join(map(str, radices)), 'signal': signal_kind,
                     'rounding': mode, 'scaling': scaling, 'dtype': dtype,
                     'preadder_round': preadder_round, 'internal_frac': internal_frac,
                     'frames': n_frames,
                     'sqnr_db': db(p_sig, p_noise),
                     'sfdr_db': (db(tone_p_sig / n_frames, tone_p_spur)
                                 if signal_kind == 'single_tone' else float('nan')),
                     'internal_peak': (max(peaks) if peaks else float('nan')),
                     'mean_err_re_lsb': bias.real, 'mean_err_im_lsb': bias.imag})
    return rows


def save_plot(rows, modes, signal_kind, dtype, scaling, out_path):
    ylabel = ('Tone SNR (dB, tone power / total error power)'
              if signal_kind == 'single_tone' else 'SQNR vs NumPy (dB)')
    fig, ax = plt.subplots(figsize=(12, 5))
    for mode in modes:
        pts = sorted((r['size'], r['sqnr_db']) for r in rows if r['rounding'] == mode)
        color, marker = MODE_STYLE.get(mode, ('#555555', 'x'))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, marker=marker,
                markersize=6, linewidth=2, label=mode)
    ax.set_xscale('log', base=2)
    sizes = sorted({r['size'] for r in rows})
    ax.set_xticks(sizes)
    ax.set_xticklabels(sizes, rotation=45, fontsize=8)
    ax.minorticks_off()
    ax.set_title(f'SQNR vs FFT size -- datapath rounding sweep ({signal_kind}, {dtype}, {scaling} scaling)')
    ax.set_xlabel('FFT size N')
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', title='rounding')
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sizes', type=int, nargs='+', default=None,
                    help=f'FFT sizes to sweep (default: quick set {QUICK_SIZES})')
    ap.add_argument('--all', action='store_true',
                    help='sweep all 53 supported top-level sizes (slow, ~an hour)')
    ap.add_argument('--signal', default='single_tone',
                    choices=['single_tone', 'sine', 'ramp', 'multitone', 'qam16'])
    ap.add_argument('--rounding', nargs='+', default=ROUNDING_MODES,
                    choices=ROUNDING_MODES, help='modes to compare')
    ap.add_argument('--dtype', default='fxp-s18/16', help='datapath dtype (RTL: fxp-s18/16)')
    ap.add_argument('--twiddle-dtype', default='fxp-s18/16', help='twiddle dtype (RTL: fxp-s18/16)')
    ap.add_argument('--coeff-dtype', default='fxp-s18/16', help='coefficient dtype (RTL: fxp-s18/16)')
    ap.add_argument('--scaling', default='fixed', choices=['fixed', 'balanced'],
                    help="per-stage shift schedule: 'fixed' = ceil(log2(radix)) always, "
                         "'balanced' = shifts chosen per config so the level rides just "
                         "under full scale (radix-3: 1-2 bits, radix-5: 2-3 bits)")
    ap.add_argument('--preadder-round', default='exit', choices=['node', 'exit'],
                    help="'exit' (default, the RTL convention) = 18 bits in memory, wide "
                         "in flight: preadder internals full precision, one shift+round "
                         "at the outputs; 'node' = legacy per-intermediate "
                         "re-quantization (pre-migration RTL)")
    ap.add_argument('--internal-frac', type=int, default=22,
                    help='exit mode: fraction bits kept on the preadder coefficient '
                         'products before the t3 adds (default 22 = c_fxp_prod_frac_width; '
                         'None-like full precision via a large value)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--frames', type=int, default=4,
                    help='frames (seeds) averaged per point for random signals '
                         '(qam16); deterministic signals always run 1 frame')
    ap.add_argument('--jobs', type=int, default=0,
                    help='parallel worker processes; 0 (default) = one per size, '
                         'capped at half the CPUs; 1 = sequential')
    ap.add_argument('--out-dir', default=str(Path(__file__).parent / 'results'))
    return ap.parse_args()


def main():
    args = parse_args()
    sizes = args.sizes or (supported_fft_sizes() if args.all else QUICK_SIZES)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    worker = partial(sweep_size, signal_kind=args.signal, modes=args.rounding,
                     dtype=args.dtype, twiddle_dtype=args.twiddle_dtype,
                     coeff_dtype=args.coeff_dtype, seed=args.seed,
                     scaling=args.scaling, preadder_round=args.preadder_round,
                     internal_frac=args.internal_frac, frames=args.frames)
    jobs = args.jobs or min(len(sizes), max(1, (os.cpu_count() or 2) // 2))

    rows = []

    def report(size_rows):
        rows.extend(size_rows)
        if size_rows[0]['signal'] == 'single_tone':
            cells = [f"{r['rounding']}: snr={r['sqnr_db']:6.2f} sfdr={r['sfdr_db']:6.2f}"
                     for r in size_rows]
        else:
            cells = [f"{r['rounding']}={r['sqnr_db']:6.2f} dB" for r in size_rows]
        print(f"N={size_rows[0]['size']:5d} ({size_rows[0]['radices']:>22}): " +
              '  '.join(cells), flush=True)

    if jobs > 1:
        with Pool(processes=jobs) as pool:
            for size_rows in pool.imap(worker, sizes):
                report(size_rows)
    else:
        for n in sizes:
            report(worker(n))

    csv_path = out_dir / 'sqnr_rounding_sweep.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    png_path = out_dir / 'sqnr_rounding_sweep.png'
    save_plot(rows, args.rounding, args.signal, args.dtype, args.scaling, png_path)

    print(f'\nwrote {csv_path}\nwrote {png_path}')


if __name__ == '__main__':
    main()
