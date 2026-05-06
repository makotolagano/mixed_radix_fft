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
    # allow a single config (int) or a per-stage iterable (list/tuple/ndarray)
    if hasattr(config, '__iter__') and not isinstance(config, (str, bytes)):
        configs = list(config)
    else:
        configs = [config] * len(stage_sizes)

    if len(configs) != len(stage_sizes):
        raise ValueError('config must be scalar or have same length as stage_sizes')

    stages_fp = []
    stages_fxp = []

    for cfg, size in zip(configs, stage_sizes):
        if cfg == 5:
            delay = size // 5
        elif cfg == 3:
            delay = size // 3
        else:
            delay = size // 2

        stages_fp.append(
            MixedRadix_SDF_stage_counter_ctrl(
                config=cfg,
                stage_index=0,
                size=size,
                cfg_delay=delay,
            )
        )
        stages_fxp.append(
            MixedRadix_SDF_stage_counter_ctrl_FXP(
                config=cfg,
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

    preadder_widths = {}
    for stage_index, stage in enumerate(stages_fxp):
        if hasattr(stage, 'pre_adder') and hasattr(stage.pre_adder, 'get_max_width_trace'):
            preadder_widths[f'stage_{stage_index}'] = stage.pre_adder.get_max_width_trace()

    return np.array(out_fp), np.array(out_fxp), preadder_widths


def sqnr(ref, test):
    noise = ref - test
    p_signal = np.mean(np.abs(ref) ** 2)
    p_noise = np.mean(np.abs(noise) ** 2)
    if p_noise == 0:
        return np.inf
    return 10 * np.log10(p_signal / p_noise)


def generate_signal(kind, n):
    if kind == 'single_tone':
        k = 20
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


def stage_latency(config, size):
    """Calculate latency in samples for a single stage.
    
    Latency = (config - 1) * (size // config) samples
    
    Args:
        config: radix (2, 3, or 5)
        size: number of samples processed per stage
    
    Returns:
        Latency in samples
    """
    cfg_delay = size // config
    return (config - 1) * cfg_delay


def total_chain_latency(configs, stage_sizes):
    """Calculate total latency for a chain of stages.
    
    Latency is cumulative: sum of per-stage latencies.
    
    Args:
        configs: list of radix values (one per stage)
        stage_sizes: list of sizes (one per stage)
    
    Returns:
        Total latency in samples
    """
    if hasattr(configs, '__iter__') and not isinstance(configs, (str, bytes)):
        cfg_list = list(configs)
    else:
        cfg_list = [configs] * len(stage_sizes)
    
    total = 0
    for cfg, size in zip(cfg_list, stage_sizes):
        total += stage_latency(cfg, size)
    return total


def evaluate_case(case_name, config, stage_sizes, stage_radices, dtypes, signal_kind):
    n = int(np.prod(stage_radices))
    x = generate_signal(signal_kind, n)
    x_pad = np.append(x, np.zeros(n))

    fp_out, _, _ = run_chain(config=config, stage_sizes=stage_sizes, input_signal=x_pad, dtype=dtypes[0])
    latency = total_chain_latency(config, stage_sizes)
    fp_ss = fp_out[latency:latency + n]
    fp_ss = fp_ss[digit_reverse(list(reversed(stage_radices)))]

    np_ref = fft(x)
    # np_ref = np_ref[digit_reverse(stage_radices)]

    dtype_results = {}
    summary_rows = []

    for dtype in dtypes:
        _, fxp_out, preadder_widths = run_chain(config=config, stage_sizes=stage_sizes, input_signal=x_pad, dtype=dtype)
        fxp_ss = fxp_out[latency:latency + n]
        fxp_ss = fxp_ss[digit_reverse(list(reversed(stage_radices)))]

        sqnr_fp_vs_fxp = sqnr(fp_ss, fxp_ss)
        sqnr_np_vs_fxp = sqnr(np_ref, fxp_ss)
        mse_fp_vs_fxp = np.mean(np.abs(fp_ss - fxp_ss) ** 2)

        dtype_results[dtype] = {
            'fft_fxp': fxp_ss,
            'sqnr_fp_vs_fxp': sqnr_fp_vs_fxp,
            'sqnr_np_vs_fxp': sqnr_np_vs_fxp,
            'mse_fp_vs_fxp': mse_fp_vs_fxp,
            'preadder_widths': preadder_widths,
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
    ax1.plot(bins, ref_mag_db, 'k-', linewidth=1.8, label='NumPy reference')
    ax1.plot(bins, fp_mag_db, 'g--', linewidth=1.8, label='FP model')

    for dtype, result in case_result['dtype_results'].items():
        fxp_mag_db = 20 * np.log10(np.maximum(np.abs(result['fft_fxp']), 1e-12))
        ax1.plot(bins, fxp_mag_db, 'r', linewidth=1.2, label=f'FXP {dtype}')

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


def save_preadder_width_plot(case_result, out_dir):
    plot_paths = []
    case_tag = case_result['case_name'].lower().replace(' ', '_').replace('(', '').replace(')', '').replace('-', '_')

    for dtype, result in case_result['dtype_results'].items():
        stage_maps = result.get('preadder_widths', {})
        if not stage_maps:
            continue

        stage_names = list(stage_maps.keys())
        operand_names = sorted({name for stage_data in stage_maps.values() for name in stage_data.keys()})

        heat_real = np.zeros((len(operand_names), len(stage_names)))
        heat_imag = np.zeros((len(operand_names), len(stage_names)))

        for col, stage_name in enumerate(stage_names):
            stage_data = stage_maps[stage_name]
            for row, operand_name in enumerate(operand_names):
                entry = stage_data.get(operand_name)
                if entry is None:
                    continue
                heat_real[row, col] = entry['real']['n_word']
                heat_imag[row, col] = entry['imag']['n_word']

        fig, axes = plt.subplots(1, 2, figsize=(16, max(8, 0.25 * len(operand_names))))

        im0 = axes[0].imshow(heat_real, aspect='auto', interpolation='nearest')
        axes[0].set_title(f'Preadder Exact n_word (Real)\n{case_result["case_name"]} | {dtype}')
        axes[0].set_xticks(range(len(stage_names)))
        axes[0].set_xticklabels(stage_names, rotation=45, ha='right')
        axes[0].set_yticks(range(len(operand_names)))
        axes[0].set_yticklabels(operand_names, fontsize=8)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label='n_word (bits)')

        im1 = axes[1].imshow(heat_imag, aspect='auto', interpolation='nearest')
        axes[1].set_title(f'Preadder Exact n_word (Imag)\n{case_result["case_name"]} | {dtype}')
        axes[1].set_xticks(range(len(stage_names)))
        axes[1].set_xticklabels(stage_names, rotation=45, ha='right')
        axes[1].set_yticks(range(len(operand_names)))
        axes[1].set_yticklabels(operand_names, fontsize=8)
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label='n_word (bits)')

        fig.tight_layout()

        dtype_tag = dtype.replace('/', '_').replace('-', '_')
        fig_path = out_dir / f'{case_tag}_preadder_exact_nword_{dtype_tag}.png'
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        plot_paths.append(fig_path)

    return plot_paths


def parse_args():
    parser = argparse.ArgumentParser(
        description='SQNR + FFT visualization for mixed-radix FP vs FXP across multiple data formats.'
    )
    parser.add_argument(
        '--n-expr',
        default='5*3*2',
        help='Radix factors for N, e.g. "5*5*5*3*2*2". Factors must be 2, 3, or 5.',
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


def _parse_radices(expr):
    parts = expr.replace('x', '*').replace('X', '*').split('*')
    radices = [int(p.strip()) for p in parts if p.strip()]
    if not radices:
        raise ValueError('n-expr must contain at least one radix')
    for r in radices:
        if r not in (2, 3, 5):
            raise ValueError('n-expr factors must be 2, 3, or 5')
    return radices


def _stage_sizes_from_radices(radices):
    n = int(np.prod(radices))
    sizes = []
    acc = 1
    for r in radices:
        sizes.append(n // acc)
        acc *= r
    return sizes


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    stage_radices = _parse_radices(args.n_expr)
    stage_sizes = _stage_sizes_from_radices(stage_radices)

    cases = [
        {
            'case_name': f"Mixed-radix chain (N={'*'.join(str(r) for r in stage_radices)})",
            'config': stage_radices,
            'stage_sizes': stage_sizes,
            'stage_radices': stage_radices,
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
        all_plot_paths.extend(save_preadder_width_plot(result, out_dir))

    print('\nSaved plots:')
    for path in all_plot_paths:
        print(f'- {path}')


if __name__ == '__main__':
    main()
