#!/usr/bin/env python3
"""Isolated golden model + test-vector generator for `twiddle_rom` (octant ROM).

The RTL stores an octant of W_N and reconstructs W_N^k with conjugate/negate/swap.
Reconstruction is exact, so the value the ROM returns for (config, k) is simply the
round-to-nearest quantization of the ideal twiddle:

        W_N^k = exp(-j 2*pi k / N)          re =  cos, im = -sin
        code  = round(value * 2**FRAC_W)    (signed, saturated to the word)

This module computes that golden result TWO ways -- directly, and by replicating the
VHDL octant fold+reconstruct -- and asserts they agree, so the emitted vectors are a
faithful reference for the RTL (independent of the ROM's internal folding).

Output: a VHDL package (CFGS + SEL_WIDTH/K_WIDTH + a TESTS array of
(config_sel, k, exp_re, exp_im)) that a self-checking testbench consumes, plus a plain
`sel k re im` text file for file-driven TBs.

Format knobs MUST match mr_fft_pkg (c_twiddle_int_width / c_twiddle_frac_width).
Change them with --int/--frac if you altered the datapath (e.g. 18-bit components).
"""
import math, argparse

# ----------------------------------------------------------------------------
# Fixed-point twiddle format (defaults match mr_fft_pkg: sfixed(1 downto -15))
# ----------------------------------------------------------------------------
class Fmt:
    def __init__(self, int_w=2, frac_w=15):
        self.int_w, self.frac_w = int_w, frac_w
        self.word_w = int_w + frac_w            # bits per component
        self.scale  = 1 << frac_w               # 2**frac_w
        self.cmin   = -(1 << (self.word_w - 1))
        self.cmax   =  (1 << (self.word_w - 1)) - 1

    def sat(self, c):
        return max(self.cmin, min(self.cmax, c))

    def rnd(self, x):
        # round half away from zero, like VHDL ieee.math_real.round
        return math.floor(x + 0.5) if x >= 0 else math.ceil(x - 0.5)

    def code(self, value):
        return self.sat(int(self.rnd(value * self.scale)))

# ----------------------------------------------------------------------------
# Golden twiddle (direct) and the VHDL octant fold+reconstruct replica
# ----------------------------------------------------------------------------
def twiddle_direct(fmt, N, k):
    ang = 2.0 * math.pi * (k % N) / N
    return fmt.code(math.cos(ang)), fmt.code(-math.sin(ang))

def f_level(N):                                  # 3=octant, 2=quarter, 1=half
    return 3 if N % 4 == 0 else (2 if N % 2 == 0 else 1)

def f_bound(N):                                  # last stored octant index
    return N // 8 if N % 4 == 0 else (N // 4 if N % 2 == 0 else (N - 1) // 2)

def octant_sample(fmt, N, r):                    # what f_oct_rom stores at index r
    ang = 2.0 * math.pi * r / N
    return fmt.code(math.cos(ang)), fmt.code(-math.sin(ang))

def twiddle_via_octant(fmt, N, k):
    """Bit-accurate replica of twiddle_rom.vhd: fold k -> r, then reconstruct."""
    lv = f_level(N)
    r, m, cj = k % N, 0, 0
    if 2 * r > N:                 r = N - r;       cj ^= 1
    if lv >= 2 and 4 * r > N:     r = N // 2 - r;  m = (m + 2) % 4; cj ^= 1
    if lv >= 3 and 8 * r > N:     r = N // 4 - r;  m = (m + (1 if cj else 3)) % 4; cj ^= 1
    a, b = octant_sample(fmt, N, r)
    if cj:  b = fmt.sat(-b)
    na, nb = fmt.sat(-a), fmt.sat(-b)
    return {0: (a, b), 1: (nb, a), 2: (na, nb), 3: (b, na)}[m]

def golden(fmt, N, k):
    direct = twiddle_direct(fmt, N, k)
    octant = twiddle_via_octant(fmt, N, k)
    assert direct == octant, f"model mismatch N={N} k={k}: direct={direct} octant={octant}"
    return direct

# ----------------------------------------------------------------------------
# Vector selection: hit every fold branch (octant/quarter/half boundaries)
# ----------------------------------------------------------------------------
def k_probes(N):
    cand = {0, 1, 2, 3, N - 1, N - 2,
            N // 8, N // 8 + 1, N // 8 - 1,
            N // 4, N // 4 + 1, N // 4 - 1,
            3 * N // 8, N // 2, N // 2 + 1, N // 2 - 1,
            5 * N // 8, 3 * N // 4, 7 * N // 8}
    # a few deterministic interior points (no RNG, reproducible)
    for j in (1, 2, 3, 5, 7):
        cand.add((j * N) // 11)
    return sorted(k for k in cand if 0 <= k < N)

def clog2(n):
    r, v = 0, 1
    while v < n:
        v *= 2; r += 1
    return max(r, 1)

# ----------------------------------------------------------------------------
# Emitters
# ----------------------------------------------------------------------------
def build_vectors(fmt, configs):
    """configs: list of (radix, size).  Returns list of (sel, k, re, im)."""
    vecs = []
    for sel, (_radix, size) in enumerate(configs):
        for k in k_probes(size):
            re, im = golden(fmt, size, k)
            vecs.append((sel, k, re, im))
    return vecs

def emit_vhdl_pkg(configs, vecs, fmt):
    sel_w = clog2(len(configs))
    k_w   = clog2(max(s for _, s in configs))
    cfg = ", ".join(f"({r},{s})" for r, s in configs)
    lines = []
    lines.append("-- AUTO-GENERATED by model/twiddle_rom_ref.py -- do not edit by hand.")
    lines.append(f"-- twiddle format: sfixed({fmt.int_w-1} downto -{fmt.frac_w})  "
                 f"({fmt.word_w} b/component)")
    lines.append("library ieee; use ieee.std_logic_1164.all;")
    lines.append("library work; use work.mr_fft_pkg.all;")
    lines.append("")
    lines.append("package twiddle_vectors_pkg is")
    lines.append(f"    constant CFGS      : t_config_arr := ( {cfg} );")
    lines.append(f"    constant SEL_WIDTH : natural := {sel_w};")
    lines.append(f"    constant K_WIDTH   : natural := {k_w};")
    lines.append("    type t_tv is record sel, k, re, im : integer; end record;")
    lines.append("    type t_tv_arr is array (natural range <>) of t_tv;")
    lines.append(f"    constant TESTS : t_tv_arr := (")
    body = [f"        ({sel:>2}, {k:>4}, {re:>7}, {im:>7})" for (sel, k, re, im) in vecs]
    lines.append(",\n".join(body))
    lines.append("    );")
    lines.append("end package;")
    return "\n".join(lines) + "\n"

def emit_txt(vecs):
    # plain "sel k re im" per line for a textio-driven TB
    return "".join(f"{sel} {k} {re} {im}\n" for (sel, k, re, im) in vecs)

# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--int",  type=int, default=2,  help="c_twiddle_int_width  (default 2)")
    ap.add_argument("--frac", type=int, default=16, help="c_twiddle_frac_width (default 16)")
    ap.add_argument("--configs", default="5:300,3:243,2:3072,2:12",
                    help="comma list of radix:size (default = tb_twiddle_rom set)")
    ap.add_argument("--vhdl", default="../rtl/tb/twiddle_vectors_pkg.vhd",
                    help="output VHDL package path")
    ap.add_argument("--txt",  default="tw_vectors/twiddle_rom_vectors.txt",
                    help="output plain-text vector path")
    a = ap.parse_args()

    fmt = Fmt(a.int, a.frac)
    configs = [tuple(int(x) for x in c.split(":")) for c in a.configs.split(",")]
    vecs = build_vectors(fmt, configs)

    import os
    for path, text in ((a.vhdl, emit_vhdl_pkg(configs, vecs, fmt)), (a.txt, emit_txt(vecs))):
        os.makedirs(os.path.dirname(path), exist_ok=True) if os.path.dirname(path) else None
        with open(path, "w") as f:
            f.write(text)
        print(f"wrote {len(vecs)} vectors -> {path}")

    # quick self-report: k=0 must be W=1 for every config
    for sel, (_r, s) in enumerate(configs):
        print(f"  cfg {sel}: N={s:>5}  W^0 = {golden(fmt, s, 0)}  W^1 = {golden(fmt, s, 1)}")
