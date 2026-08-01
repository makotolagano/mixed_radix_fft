import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.fft import fft

from mixed_radix_fft import MixedRadix_SDF_stage_counter_ctrl
from mixed_radix_fft_fxp import MixedRadix_FinalScaler_FXP, MixedRadix_SDF_stage_counter_ctrl_FXP
from utils import digit_reverse

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception as exc:
    raise ImportError('matplotlib is required for plotting. Install it in your env.') from exc


def shift_schedule(radices, scaling='fixed'):
    """Per-stage preadder output right-shifts (the scaling schedule).

    'fixed'    -- every stage shifts by ceil(log2(radix)) (2->1, 3->2, 5->3).
                  Radix-3/5 stages over-scale by 4/3 and 8/5, so the level
                  decays by 2**total_shift/N through the chain and the final
                  scaler amplifies late-stage rounding noise by that factor.
    'balanced' -- radix-2 shifts 1 (exact), radix-3 shifts 1 or 2, radix-5
                  shifts 2 or 3, chosen greedily so the cumulative net gain
                  prod(radix)/prod(2**shift) never exceeds 1: the level rides
                  just under full scale the whole way and the end residual
                  stays below one bit.

    Pure INTEGER rule (prod_r <= prod_s << c) -- the VHDL elaboration-time
    schedule function must reproduce this verbatim so RTL and model derive
    identical schedules.
    """
    cands = {2: (1,), 3: (1, 2), 5: (2, 3)}
    fixed = {2: 1, 3: 2, 5: 3}
    prod_r, prod_s = 1, 1
    shifts = []
    for r in radices:
        if scaling == 'fixed':
            s = fixed[r]
        elif scaling == 'balanced':
            s = next((c for c in cands[r] if prod_r * r <= prod_s << c), cands[r][-1])
        else:
            raise ValueError(f'unknown scaling schedule: {scaling}')
        prod_r *= r
        prod_s <<= s
        shifts.append(s)
    return shifts


def run_chain(config, stage_sizes, input_signal, dtype='fxp-s32/12', twiddle_dtype=None,
              input_dtype=None, coeff_dtype=None, output_dtype=None, rounding=None,
              scaling='fixed', preadder_exit_round=False, preadder_internal_frac=None,
              run_fp=True, final_scaler=True):
    # final_scaler=False returns the RAW chain output, X * 2**(-total_shift) --
    # what the HARDWARE delivers (no output rescaler exists in the RTL); the
    # caller must then reference against fft(x) * 2**(-total_shift), NOT
    # fft(x)/N. The scaler's requantization shifts single-frame tone SQNR by
    # up to ~2.5 dB either way, so hardware-representative runs need False.
    # run_fp=False skips the floating-point reference chain entirely (the
    # sweeps compare against scipy's FFT and discard fp_out) -- ~2x faster.
    # dtype is the inner_type (datapath). input_dtype quantizes the chain input (first
    # stage only), coeff_dtype the preadder constants, output_dtype the final result.
    # rounding, when set, overrides the datapath/coefficient quantizers ('floor'
    # when unset). Twiddle table values are design-time ROM constants, so they
    # ALWAYS keep round-to-nearest regardless of the override -- matching the
    # RTL, where the ROM init is precomputed and only the datapath style would
    # ever change.
    data_rounding = 'floor' if rounding is None else rounding
    tw_rounding = 'around'
    # allow a single config (int) or a per-stage iterable (list/tuple/ndarray)
    if hasattr(config, '__iter__') and not isinstance(config, (str, bytes)):
        configs = list(config)
    else:
        configs = [config] * len(stage_sizes)

    if len(configs) != len(stage_sizes):
        raise ValueError('config must be scalar or have same length as stage_sizes')

    stages_fp = []
    stages_fxp = []
    total_shift = 0

    shifts = shift_schedule(configs, scaling)

    for stage_pos, (cfg, size) in enumerate(zip(configs, stage_sizes)):
        delay = size // cfg
        shift = shifts[stage_pos]
        total_shift += shift

        if run_fp:
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
                preadder_shift_bits=shift,
                dtype=dtype,
                twiddle_dtype=twiddle_dtype,
                coeff_dtype=coeff_dtype,
                input_dtype=(input_dtype if stage_pos == 0 else None),
                rounding=data_rounding,
                twiddle_rounding=tw_rounding,
                capability=cfg,
                preadder_exit_round=preadder_exit_round,
                preadder_internal_frac=preadder_internal_frac,
            )
        )

    out_fp = []
    out_fxp = []
    # output word defaults to the datapath word -- without this, dtype=None
    # fell through to fxpmath's DEFAULT fxp-s16/15 (range +-1, saturate),
    # silently clipping any output bin above 1.0 (caught with a sqrt(2) tone)
    out_dtype = output_dtype or dtype
    # hardware-faithful scaler: quantized s4.21 constant, wrap (mr_fft_scaler)
    scaler = MixedRadix_FinalScaler_FXP(size=stage_sizes[0], total_shift=total_shift,
                                        dtype=out_dtype, output_dtype=out_dtype,
                                        rounding=data_rounding, overflow='wrap',
                                        quantized_scale_frac=21)
    print(f"Scale = {str(2**total_shift/stage_sizes[0])}")

    for sample in input_signal:
        val_fp = sample
        val_fxp = sample

        if run_fp:
            for stage in stages_fp:
                val_fp = stage.calculate(val_fp, valid=True)
            out_fp.append(val_fp / stage_sizes[0])

        for stage in stages_fxp:
            val_fxp = stage.calculate(val_fxp, valid=True)

        if final_scaler:
            val_fxp = scaler.scale_sample(val_fxp)
        out_fxp.append(val_fxp)

    preadder_widths = {}
    for stage_index, stage in enumerate(stages_fxp):
        if hasattr(stage, 'pre_adder') and hasattr(stage.pre_adder, 'get_max_width_trace'):
            preadder_widths[f'stage_{stage_index}'] = stage.pre_adder.get_max_width_trace()
        if preadder_exit_round:
            # largest wide-datapath magnitude seen inside the preadder -- sizes
            # the internal guard bits (must fit the DSP 25-bit data port)
            preadder_widths[f'stage_{stage_index}_peak'] = stage.pre_adder.internal_peak

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
        # coherent (bin-centered) complex tone: integer k -> zero leakage, so
        # everything outside bin k is quantization noise. k ~ N/7 and coprime
        # with N so the twiddle walk is generic (a k sharing a factor with N
        # exercises only degenerate exponent subsets in some stages) -- the
        # IEEE-1241-style "relatively prime number of cycles" rule.
        # 0.9 amplitude ~ -1 dBFS, honoring the chain input contract.
        k = max(1, n // 7)
        while math.gcd(k, n) != 1:
            k += 1
        idx = np.arange(n, dtype=float)
        return np.sqrt(2) * np.exp(2j * np.pi * k * idx / n)
    if kind == 'sine':
        k = max(1, n // 7)
        idx = np.arange(n, dtype=float)
        return np.sqrt(2) * np.sin(2 * np.pi * k * idx / n)
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
    if kind == 'qam16':
        # Unit-power 16-QAM: I,Q in {-3,-1,+1,+3} scaled by 1/sqrt(10)
        const_points = np.array([-3, -1, 1, 3], dtype=float) / np.sqrt(10)
        i_idx = np.random.randint(0, 4, size=n)
        q_idx = np.random.randint(0, 4, size=n)
        return const_points[i_idx] + 1j * const_points[q_idx]
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


def evaluate_case(case_name, config, stage_sizes, stage_radices, dtypes, twiddle_dtypes, signal_kind,
                  input_dtype=None, coeff_dtype=None, output_dtype=None, rounding=None,
                  scaling='fixed', preadder_round='exit', internal_frac=22):
    n = int(np.prod(stage_radices))
    x = generate_signal(signal_kind, n)
    x_pad = np.append(x, np.zeros(n))

    fp_out, _, _ = run_chain(config=config, stage_sizes=stage_sizes, input_signal=x_pad, dtype=dtypes[0], twiddle_dtype=(twiddle_dtypes[0] if twiddle_dtypes else None))
    latency = total_chain_latency(config, stage_sizes)
    fp_ss = fp_out[latency:latency + n]
    fp_ss = fp_ss[digit_reverse(list(reversed(stage_radices)))]

    np_ref = fft(x) / n
    # np_ref = np_ref[digit_reverse(stage_radices)]

    dtype_results = {}
    summary_rows = []

    if not twiddle_dtypes:
        twiddle_dtypes = [None]

    for dtype in dtypes:
        for twiddle_dtype in twiddle_dtypes:
            _, fxp_out, preadder_widths = run_chain(
                config=config,
                stage_sizes=stage_sizes,
                input_signal=x_pad,
                dtype=dtype,
                twiddle_dtype=twiddle_dtype,
                input_dtype=input_dtype,
                coeff_dtype=coeff_dtype,
                output_dtype=output_dtype,
                rounding=rounding,
                scaling=scaling,
                preadder_exit_round=(preadder_round == 'exit'),
                preadder_internal_frac=internal_frac,
                run_fp=False,   # the float reference was already computed above
            )
            fxp_ss = fxp_out[latency:latency + n]
            fxp_ss = fxp_ss[digit_reverse(list(reversed(stage_radices)))]

            sqnr_fp_vs_fxp = sqnr(fp_ss, fxp_ss)
            sqnr_np_vs_fxp = sqnr(np_ref, fxp_ss)
            mse_np_vs_fxp = np.mean(np.abs(np_ref - fxp_ss) ** 2)

            label = f'{dtype} | twiddle={twiddle_dtype or dtype}'
            dtype_results[label] = {
                'fft_fxp': fxp_ss,
                'sqnr_fp_vs_fxp': sqnr_fp_vs_fxp,
                'sqnr_np_vs_fxp': sqnr_np_vs_fxp,
                'mse_np_vs_fxp': mse_np_vs_fxp,
                'preadder_widths': preadder_widths,
            }
            summary_rows.append((label, sqnr_fp_vs_fxp, sqnr_np_vs_fxp, mse_np_vs_fxp))

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
        # keep only the per-operand width-trace dicts; the flat 'stage_N_peak'
        # floats (exit-mode internal peak, for guard-bit sizing) are not part
        # of the heatmap
        stage_maps = {k: v for k, v in result.get('preadder_widths', {}).items()
                      if isinstance(v, dict)}
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
        default=['fxp-s18/16'],
        help='List of fxpmath dtype formats to compare (default: the RTL data word).',
    )
    parser.add_argument(
        '--signal',
        choices=['single_tone', 'ramp', 'sine', 'multitone', 'qam16'],
        default='qam16',
        help='Input signal used for FFT test.',
    )
    parser.add_argument(
        '--out-dir',
        default='results',
        help='Directory for saved plot images.',
    )
    parser.add_argument(
        '--twiddle_dtype',
        default=None,
        help='Twiddle factor dtype (e.g. "fxp-s32/24", "fxp-s64/48"). Defaults to --dtypes if not specified.',
    )
    parser.add_argument(
        '--input_dtype',
        default=None,
        help='Chain input word (fxpmath dtype). Defaults to the datapath --dtypes value.',
    )
    parser.add_argument(
        '--coeff_dtype',
        default=None,
        help='Preadder coefficient word. Defaults to the datapath --dtypes value.',
    )
    parser.add_argument(
        '--output_dtype',
        default=None,
        help='Final result word out of the scaler. Defaults to the datapath --dtypes value.',
    )
    parser.add_argument(
        '--rounding',
        choices=['floor', 'around', 'half_up'],
        default='half_up',
        help="Rounding for the datapath and coefficient quantizers (default 'half_up', "
             "the RTL convention). 'around' = ties-to-even (VHDL fixed_round), "
             "'floor' = truncation (legacy). Twiddle table values are ROM constants "
             "and always stay round-to-nearest.",
    )
    parser.add_argument(
        '--scaling',
        choices=['fixed', 'balanced'],
        default='fixed',
        help="Per-stage shift schedule (default 'fixed' = the RTL: ceil(log2(radix)) "
             "derived from the radix mode; 'balanced' is model-only for now).",
    )
    parser.add_argument(
        '--preadder-round',
        choices=['node', 'exit'],
        default='exit',
        help="'exit' (default, the RTL convention) = 18 bits in memory, wide in "
             "flight, one shift+round at the preadder outputs; 'node' = legacy "
             "per-intermediate re-quantization.",
    )
    parser.add_argument(
        '--internal-frac',
        type=int,
        default=22,
        help='exit mode: fraction bits kept on the preadder products '
             '(c_fxp_prod_frac_width; default 22, the RTL value).',
    )
    parser.add_argument(
        '--sweep-file',
        default='sweep_config.json',
        help='JSON file containing sweep settings for dtypes and twiddle_dtypes.',
    )
    return parser.parse_args()


def _load_sweep_config(path):
    sweep_path = Path(path)
    if not sweep_path.exists():
        return {}
    with sweep_path.open('r', encoding='utf-8') as f:
        return json.load(f)


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

    sweep_cfg = _load_sweep_config(args.sweep_file)
    if sweep_cfg:
        args.n_expr = sweep_cfg.get('n_expr', args.n_expr)
        args.dtypes = sweep_cfg.get('dtypes', args.dtypes)
        args.signals = sweep_cfg.get('signals', [sweep_cfg.get('signal', args.signal)])
        if 'twiddle_dtypes' in sweep_cfg:
            args.twiddle_dtypes = sweep_cfg.get('twiddle_dtypes')
        else:
            args.twiddle_dtypes = [args.twiddle_dtype] if args.twiddle_dtype else [None]
        args.input_dtype = sweep_cfg.get('input_dtype', args.input_dtype)
        args.coeff_dtype = sweep_cfg.get('coeff_dtype', args.coeff_dtype)
        args.output_dtype = sweep_cfg.get('output_dtype', args.output_dtype)
        args.rounding = sweep_cfg.get('rounding', args.rounding)
        args.scaling = sweep_cfg.get('scaling', args.scaling)
        args.preadder_round = sweep_cfg.get('preadder_round', args.preadder_round)
        args.internal_frac = sweep_cfg.get('internal_frac', args.internal_frac)
    else:
        args.signals = [args.signal]
        args.twiddle_dtypes = [args.twiddle_dtype] if args.twiddle_dtype else [None]

    stage_radices = _parse_radices(args.n_expr)
    stage_sizes = _stage_sizes_from_radices(stage_radices)

    cases = [
        {
            'case_name': f"Mixed-radix chain (N={'*'.join(str(r) for r in stage_radices)}) | signal={signal_kind}",
            'config': stage_radices,
            'stage_sizes': stage_sizes,
            'stage_radices': stage_radices,
            'signal_kind': signal_kind,
        }
        for signal_kind in args.signals
    ]

    all_plot_paths = []
    for case in cases:
        result = evaluate_case(
            case_name=case['case_name'],
            config=case['config'],
            stage_sizes=case['stage_sizes'],
            stage_radices=case['stage_radices'],
            dtypes=args.dtypes,
            twiddle_dtypes=args.twiddle_dtypes,
            signal_kind=case['signal_kind'],
            input_dtype=args.input_dtype,
            coeff_dtype=args.coeff_dtype,
            output_dtype=args.output_dtype,
            rounding=args.rounding,
            scaling=args.scaling,
            preadder_round=args.preadder_round,
            internal_frac=args.internal_frac,
        )
        print_summary_table(result)
        all_plot_paths.extend(save_plots(result, out_dir))
        all_plot_paths.extend(save_preadder_width_plot(result, out_dir))

    print('\nSaved plots:')
    for path in all_plot_paths:
        print(f'- {path}')


if __name__ == '__main__':
    main()
