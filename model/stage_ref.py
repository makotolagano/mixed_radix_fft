"""Golden-vector generator for rtl/src/mr_fft_stage.vhd.

Runs MixedRadix_SDF_stage_counter_ctrl_FXP with the RTL's fixed-point
conventions: "18 bits in memory, wide in flight" -- data s18/16, coeff/
twiddle s18/16, wrap overflow, half_up rounding, preadder exit_round with
internal_frac=22 (docs/datapath_width_convention.md). --cap selects the
hardware variant (RTL G_CAPABILITY 2/1/0 = model capability 5/3/2):

  * octant-ROM twiddles (TwiddleOctantROM == rtl/src/twiddle_rom.vhd, verified
    bit-exact by model/twiddle_rom_ref.py + tb_twiddle_rom_gen): output beat
    (arm p, position k) is rotated by W_N^(p*k*radix^STAGE_INDEX)
  * STAGE_INDEX = 0 (first stage of a chain)
  * preadder_shift_bits = None -- the FIXED scaling schedule, i.e. the
    per-radix default ceil(log2(radix)) that the RTL preadder derives from
    its s0/s1 mode and fuses into the half-up exit rounding

The runtime delay (size/radix) is a stage port (i_config_delay), so configs
with different delays are exercised: three radices at delay 4 plus a radix-2
config at delay 2 (< the physical FIFO depth, which is the max delay over
G_CONFIGS).

Per config it emits the full input stream (2 frames + 1 zero flush frame) and
the expected output stream (valid from call index N - N/radix) as raw
two's-complement integer codes.

Regenerate:  .venv/bin/python model/stage_ref.py --cap 2   (and --cap 1, --cap 0)
Consumed by: rtl/tb/tb_mr_fft_stage.vhd (cap2),
             rtl/tb/tb_mr_fft_stage_cap1.vhd, rtl/tb/tb_mr_fft_stage_cap0.vhd
"""
import argparse

import numpy as np

from mixed_radix_fft_fxp import MixedRadix_SDF_stage_counter_ctrl_FXP

# per RTL G_CAPABILITY: (radix, N) configs; keep configs with a smaller delay
# than the max (runtime delay < physical FIFO depth) and delay-1 configs
# (exercise the FIFO depth-1 bypass register).
# Deliberately shuffled delay order (grow AND shrink transitions): the TBs
# reconfigure on the fly via the i_reconfig strobe, which must make the order
# irrelevant.
CAP_CONFIGS = {
    2: [(2, 8), (2, 2), (5, 20), (2, 4), (3, 12), (5, 5)],  # delays 4,1,4,2,4,1
    1: [(2, 8), (3, 3), (3, 12), (2, 4), (3, 6), (2, 2)],   # delays 4,1,4,2,2,1
    0: [(2, 8), (2, 2), (2, 4)],                            # delays 4,1,2
}
MODEL_CAP = {2: 5, 1: 3, 0: 2}               # RTL G_CAPABILITY -> model capability
STAGE_INDEX = 0                              # RTL G_STAGE_INDEX
TWIDDLE_DTYPE = 'fxp-s18/16'                 # c_twiddle_int/frac_width
N_FRAMES = 2


def to_code(value, frac_w):
    code = round(value * 2.0 ** frac_w)
    assert code == value * 2.0 ** frac_w, f'non-exact code for {value}'
    return int(code)


def input_codes(radix, size, word_w, rng):
    """2 frames: a small deterministic ramp, then seeded random (kept within
    +-2**(word_w-4) so radix-5 sums stay inside the data word -- wrap behavior
    is already covered by the preadder TB), then a zero flush frame."""
    ramp = [(64 * (k + 1), -32 * (k + 1)) for k in range(size)]
    lim = 2 ** (word_w - 4)
    rand = [(int(rng.integers(-lim, lim)), int(rng.integers(-lim, lim)))
            for _ in range(size)]
    return ramp + rand + [(0, 0)] * size


def run_stage(radix, size, codes, data_dtype, coeff_dtype, frac_w, capability):
    stage = MixedRadix_SDF_stage_counter_ctrl_FXP(
        config=radix, stage_index=STAGE_INDEX, size=size, cfg_delay=size // radix,
        dtype=data_dtype, coeff_dtype=coeff_dtype, preadder_shift_bits=None,
        overflow='wrap', rounding='half_up', capability=capability,
        twiddle_source='octant', twiddle_dtype=TWIDDLE_DTYPE,
        preadder_exit_round=True, preadder_internal_frac=22)
    scale = 2.0 ** frac_w
    out = []
    for (re, im) in codes:
        y = stage.calculate(complex(re / scale, im / scale))
        out.append((to_code(y.real, frac_w), to_code(y.imag, frac_w)))
    return out


def build(seed, word_w, data_dtype, coeff_dtype, frac_w, cap):
    rng = np.random.default_rng(seed)
    cfgs = []
    for (radix, size) in CAP_CONFIGS[cap]:
        delay = size // radix
        cin = input_codes(radix, size, word_w, rng)
        cout = run_stage(radix, size, cin, data_dtype, coeff_dtype, frac_w,
                         MODEL_CAP[cap])
        valid_start = size - delay
        exp = cout[valid_start:valid_start + N_FRAMES * size]
        cfgs.append(dict(radix=radix, size=size, valid_start=valid_start,
                         cin=cin, exp=exp))
    return cfgs


def fmt_arr(pairs):
    """Interleave (re, im) pairs into a VHDL integer aggregate, 8 per line."""
    flat = [v for p in pairs for v in p]
    lines = []
    for i in range(0, len(flat), 8):
        lines.append(", ".join(f"{v:>7}" for v in flat[i:i + 8]))
    return ",\n        ".join(lines)


def emit_vhdl_pkg(cfgs, word_w, frac_w, pkg_name, cap):
    in_ofs, exp_ofs = [0], [0]
    for c in cfgs:
        in_ofs.append(in_ofs[-1] + 2 * len(c['cin']))
        exp_ofs.append(exp_ofs[-1] + 2 * len(c['exp']))

    L = []
    L.append(f"-- AUTO-GENERATED by model/stage_ref.py --cap {cap} -- do not edit by hand.")
    L.append(f"-- data format: sfixed({word_w - frac_w - 1} downto -{frac_w})"
             f"  ({word_w} b/component), wrap + round-to-nearest")
    L.append("-- Streams are interleaved integer codes: re0, im0, re1, im1, ...")
    L.append(f"package {pkg_name} is")
    L.append("    type t_int_arr is array (natural range <>) of integer;")
    L.append(f"    constant NUM_CFGS : integer := {len(cfgs)};")
    L.append("    constant CFG_RADIX : t_int_arr := ("
             + ", ".join(str(c['radix']) for c in cfgs) + ");")
    L.append("    constant CFG_SIZE : t_int_arr := ("
             + ", ".join(str(c['size']) for c in cfgs) + ");")
    L.append("    -- model call index of the first valid output (= size - delay)")
    L.append("    constant CFG_VALID_START : t_int_arr := ("
             + ", ".join(str(c['valid_start']) for c in cfgs) + ");")
    L.append("    -- slices of IN_CODES / EXP_CODES per config: [ofs(i), ofs(i+1))")
    L.append("    constant CFG_IN_OFS : t_int_arr := ("
             + ", ".join(str(v) for v in in_ofs) + ");")
    L.append("    constant CFG_EXP_OFS : t_int_arr := ("
             + ", ".join(str(v) for v in exp_ofs) + ");")
    L.append("    constant IN_CODES : t_int_arr := (")
    L.append("        " + ",\n        ".join(fmt_arr(c['cin']) for c in cfgs))
    L.append("    );")
    L.append("    constant EXP_CODES : t_int_arr := (")
    L.append("        " + ",\n        ".join(fmt_arr(c['exp']) for c in cfgs))
    L.append("    );")
    L.append("end package;")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--int", type=int, default=2, help="c_fxp_int_width (default 2)")
    ap.add_argument("--frac", type=int, default=16, help="c_fxp_frac_width (default 16)")
    ap.add_argument("--coeff-int", type=int, default=2)
    ap.add_argument("--coeff-frac", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cap", type=int, default=2, choices=[2, 1, 0],
                    help="RTL G_CAPABILITY (default 2 = radix235)")
    ap.add_argument("--vhdl", default=None,
                    help="output file (default rtl/tb/stage_vectors[_capN]_pkg.vhd)")
    a = ap.parse_args()

    # cap 2 keeps the original package name/path for compatibility
    pkg_name = "stage_vectors_pkg" if a.cap == 2 else f"stage_vectors_cap{a.cap}_pkg"
    if a.vhdl is None:
        a.vhdl = f"rtl/tb/{pkg_name}.vhd"

    word_w = a.int + a.frac
    data_dtype = f"fxp-s{word_w}/{a.frac}"
    coeff_dtype = f"fxp-s{a.coeff_int + a.coeff_frac}/{a.coeff_frac}"

    cfgs = build(a.seed, word_w, data_dtype, coeff_dtype, a.frac, a.cap)
    with open(a.vhdl, "w") as f:
        f.write(emit_vhdl_pkg(cfgs, word_w, a.frac, pkg_name, a.cap))
    total_in = sum(len(c['cin']) for c in cfgs)
    total_exp = sum(len(c['exp']) for c in cfgs)
    print(f"wrote {len(cfgs)} configs ({total_in} input samples, "
          f"{total_exp} expected samples) -> {a.vhdl}")
    for c in cfgs:
        print(f"  radix={c['radix']} N={c['size']} valid_start={c['valid_start']}"
              f" n_in={len(c['cin'])} n_exp={len(c['exp'])}")
