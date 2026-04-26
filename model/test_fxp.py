import argparse
import numpy as np
from scipy.fft import fft

from mixed_radix_fft import MixedRadix_SDF_stage_counter_ctrl
from mixed_radix_fft_fxp import MixedRadix_SDF_stage_counter_ctrl_FXP
from utils import digit_reverse


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


def _fmt_complex(val):
    return f"{val.real:>10.4f} {val.imag:+10.4f}j"


def print_result_table(ss_fp, ss_fxp, ref_fft, rows=8):
    rows = min(rows, len(ss_fp), len(ss_fxp), len(ref_fft))
    line = "-" * 95
    print(line)
    print(f"{'idx':>4} | {'FP output':>25} | {'FXP output':>25} | {'NumPy FFT':>25}")
    print(line)
    for idx in range(rows):
        print(
            f"{idx:>4} | "
            f"{_fmt_complex(ss_fp[idx]):>25} | "
            f"{_fmt_complex(ss_fxp[idx]):>25} | "
            f"{_fmt_complex(ref_fft[idx]):>25}"
        )
    print(line)


def report_case(name, config, stage_sizes, stage_radices, N, dtype):
    x = np.arange(N).astype(float)
    x_pad = np.append(x, np.zeros(N))

    out_fp, out_fxp = run_chain(
        config=config,
        stage_sizes=stage_sizes,
        input_signal=x_pad,
        dtype=dtype,
    )

    ss_fp = out_fp[N - 1:N - 1 + N]
    ss_fxp = out_fxp[N - 1:N - 1 + N]

    ref_fft = fft(x)
    ref_fft = ref_fft[digit_reverse(stage_radices)]

    print(f'\n=== {name} ===')
    print(f'Config={config} | stage_sizes={stage_sizes} | stage_radices={stage_radices} | N={N} | dtype={dtype}')
    print_result_table(ss_fp, ss_fxp, ref_fft, rows=8)
    print(f'SQNR(FP vs FXP): {sqnr(ss_fp, ss_fxp):.2f} dB')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mixed-radix FXP validation cases.')
    parser.add_argument(
        '--dtype',
        default='fxp-s32/12',
        help="Fixed-point dtype passed to fxpmath (e.g. 'fxp-s32/12', 'fxp-s24/10').",
    )
    args = parser.parse_args()

    report_case(
        name='Radix-3 chain (9 = 3x3)',
        config=3,
        stage_sizes=[27, 9, 3],
        stage_radices=[3, 3, 3],
        N=27,
        dtype=args.dtype,
    )
