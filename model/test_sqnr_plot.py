import argparse
from pathlib import Path

import numpy as np
from scipy.fft import fft

from mixed_radix_fft import MixedRadix_SDF_stage_counter_ctrl
from mixed_radix_fft_fxp import MixedRadix_SDF_stage_counter_ctrl_FXP
from utils import digit_reverse

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception as exc:
    raise ImportError('matplotlib is required for plotting. Install it in your env.') from exc


def run_chain(config, stage_sizes, input_signal, dtype='fxp-s32/12'):
    stages_fp = []
    stages_fxp = []

    for size in stage_sizes:
        if config == 5:
            delay = size // 5
        elif config == 3:
            delay = size // 3
        else:
            delay = size // 2

        stages_fp.append(
            MixedRadix_SDF_stage_counter_ctrl(
                config=config,
                stage_index=0,
                size=size,
                cfg_delay=delay,
            )
        )
        stages_fxp.append(
            MixedRadix_SDF_stage_counter_ctrl_FXP(
                config=config,
                stage_index=0,
                size=size,
                cfg_delay=delay,
                dtype=dtype,
            )
        )

    out_fp = []
    out_fxp = []

    for sample in input_signal:
        val_fp = sample
        val_fxp = sample

        for stage in stages_fp:
            val_fp = stage.calculate(val_fp, valid=True)

        for stage in stages_fxp:
            val_fxp = stage.calculate(val_fxp, valid=True)

        out_fp.append(val_fp)
        out_fxp.append(val_fxp)

    return np.array(out_fp), np.array(out_fxp)


def sqnr(ref, test):
    noise = ref - test
    p_signal = np.mean(np.abs(ref) ** 2)
    p_noise = np.mean(np.abs(noise) ** 2)
    if p_noise == 0:
        return np.inf
    return 10 * np.log10(p_signal / p_noise)


def generate_signal(kind, n):
    if kind == 'single_tone':
        k = 7
        idx = np.arange(n, dtype=float)
        return 0.9 * np.sin(2 * np.pi * k * idx / n)
    if kind == 'sine':
        k = max(1, n // 7)
        idx = np.arange(n, dtype=float)
        return 0.9 * np.sin(2 * np.pi * k * idx / n)
    if kind == 'ramp':
        return np.arange(n, dtype=float)
    if kind == 'multitone':
        idx = np.arange(n, dtype=float)
        tones = [1, max(2, n // 8), max(3, n // 5)]
        sig = np.zeros(n, dtype=float)
        for tone in tones:
            sig += np.sin(2 * np.pi * tone * idx / n)
        sig = 0.9 * sig / np.max(np.abs(sig))
        return sig
    raise ValueError(f'Unsupported signal kind: {kind}')


def evaluate_case(case_name, config, stage_sizes, stage_radices, dtypes, signal_kind):
    n = int(np.prod(stage_radices))
    x = generate_signal(signal_kind, n)
    x_pad = np.append(x, np.zeros(n))

    fp_out, _ = run_chain(config=config, stage_sizes=stage_sizes, input_signal=x_pad, dtype=dtypes[0])
    fp_ss = fp_out[n - 1:n - 1 + n]

    np_ref = fft(x)
    np_ref = np_ref[digit_reverse(stage_radices)]

    dtype_results = {}
    summary_rows = []

    for dtype in dtypes:
        _, fxp_out = run_chain(config=config, stage_sizes=stage_sizes, input_signal=x_pad, dtype=dtype)
        fxp_ss = fxp_out[n - 1:n - 1 + n]

        sqnr_fp_vs_fxp = sqnr(fp_ss, fxp_ss)
        sqnr_np_vs_fxp = sqnr(np_ref, fxp_ss)
        mse_fp_vs_fxp = np.mean(np.abs(fp_ss - fxp_ss) ** 2)

        dtype_results[dtype] = {
            'fft_fxp': fxp_ss,
            'sqnr_fp_vs_fxp': sqnr_fp_vs_fxp,
            'sqnr_np_vs_fxp': sqnr_np_vs_fxp,
            'mse_fp_vs_fxp': mse_fp_vs_fxp,
        }
        summary_rows.append((dtype, sqnr_fp_vs_fxp, sqnr_np_vs_fxp, mse_fp_vs_fxp))

    return {
        'case_name': case_name,
        'n': n,
        'signal': x,
        'np_ref': np_ref,
        'fp_ss': fp_ss,
        'dtype_results': dtype_results,
        'summary_rows': summary_rows,
    }


def print_summary_table(case_result):
    print(f"\n=== {case_result['case_name']} (N={case_result['n']}) ===")
    print('-' * 88)
    print(f"{'dtype':<18} | {'SQNR FP-FXP (dB)':>16} | {'SQNR NP-FXP (dB)':>16} | {'MSE FP-FXP':>16}")
    print('-' * 88)
    for dtype, sqnr_fp_fxp, sqnr_np_fxp, mse in case_result['summary_rows']:
        print(f"{dtype:<18} | {sqnr_fp_fxp:>16.4f} | {sqnr_np_fxp:>16.4f} | {mse:>16.6e}")
    print('-' * 88)


def save_plots(case_result, out_dir):
    case_tag = case_result['case_name'].lower().replace(' ', '_').replace('(', '').replace(')', '').replace('-', '_')
    bins = np.arange(case_result['n'])

    ref_mag_db = 20 * np.log10(np.maximum(np.abs(case_result['np_ref']), 1e-12))
    fp_mag_db = 20 * np.log10(np.maximum(np.abs(case_result['fp_ss']), 1e-12))

    fig1, ax1 = plt.subplots(figsize=(12, 5))
    ax1.plot(bins, ref_mag_db, 'k-', linewidth=2.0, label='NumPy reference')
    ax1.plot(bins, fp_mag_db, 'b--', linewidth=1.8, label='FP model')

    for dtype, result in case_result['dtype_results'].items():
        fxp_mag_db = 20 * np.log10(np.maximum(np.abs(result['fft_fxp']), 1e-12))
        ax1.plot(bins, fxp_mag_db, linewidth=1.2, label=f'FXP {dtype}')

    ax1.set_title(f"FFT Magnitude Comparison - {case_result['case_name']}")
    ax1.set_xlabel('FFT bin')
    ax1.set_ylabel('Magnitude (dB)')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=8)
    fig1.tight_layout()
    fig1_path = out_dir / f'{case_tag}_fft_magnitude.png'
    fig1.savefig(fig1_path, dpi=160)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(12, 5))
    for dtype, result in case_result['dtype_results'].items():
        err_db = 20 * np.log10(np.maximum(np.abs(case_result['fp_ss'] - result['fft_fxp']), 1e-12))
        ax2.plot(bins, err_db, linewidth=1.2, label=f'|FP-FXP| {dtype}')

    ax2.set_title(f"Per-Bin Error vs FP - {case_result['case_name']}")
    ax2.set_xlabel('FFT bin')
    ax2.set_ylabel('Error magnitude (dB)')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=8)
    fig2.tight_layout()
    fig2_path = out_dir / f'{case_tag}_error_vs_fp.png'
    fig2.savefig(fig2_path, dpi=160)
    plt.close(fig2)

    return [fig1_path, fig2_path]


def parse_args():
    parser = argparse.ArgumentParser(
        description='SQNR + FFT visualization for mixed-radix FP vs FXP across multiple data formats.'
    )
    parser.add_argument(
        '--dtypes',
        nargs='+',
        default=['fxp-s16/10', 'fxp-s24/10', 'fxp-s32/12'],
        help='List of fxpmath dtype formats to compare.',
    )
    parser.add_argument(
        '--signal',
        choices=['single_tone', 'ramp', 'sine', 'multitone'],
        default='single_tone',
        help='Input signal used for FFT test.',
    )
    parser.add_argument(
        '--out-dir',
        default='results',
        help='Directory for saved plot images.',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        {
            'case_name': 'Radix-3 chain (243 = 3^5)',
            'config': 3,
            'stage_sizes': [243, 81, 27, 9, 3],
            'stage_radices': [3, 3, 3, 3, 3],
        }
    ]

    all_plot_paths = []
    for case in cases:
        result = evaluate_case(
            case_name=case['case_name'],
            config=case['config'],
            stage_sizes=case['stage_sizes'],
            stage_radices=case['stage_radices'],
            dtypes=args.dtypes,
            signal_kind=args.signal,
        )
        print_summary_table(result)
        all_plot_paths.extend(save_plots(result, out_dir))

    print('\nSaved plots:')
    for path in all_plot_paths:
        print(f'- {path}')


if __name__ == '__main__':
    main()
