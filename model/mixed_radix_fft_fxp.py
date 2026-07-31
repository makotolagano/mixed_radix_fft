import math
import re
from collections import deque

import numpy as np
from fxpmath import Fxp
from utils import Fifo

# When True, every quantization goes through fxpmath (the original reference
# implementation). The default fast path does the same arithmetic in plain
# float64/int -- bit-identical (all reachable values are exact in float64;
# proven by the golden-vector regeneration and model/check_fast_quantizer) --
# and ~50-100x faster, since Fxp() deep-copies its config on every call.
USE_REFERENCE_FXP = False


class _Quantizer:
    def __init__(self, dtype='fxp-s32/12', overflow='wrap', rounding='floor'):
        # rounding='floor' models VHDL fixed_truncate (toward -inf); 'around' is
        # round-to-nearest ties-to-even (VHDL fixed_round). 'half_up' is
        # round-to-nearest ties-toward-+inf: +LSB/2 then floor, mirroring the
        # DSP-post-adder-friendly RTL form `resize(x + half_lsb, ...,
        # fixed_truncate)`. Not a native fxpmath mode, so it is applied as a
        # pre-bias in q()/qw(). Chosen per-quantizer so it can be set globally.
        self._half_up = (rounding == 'half_up')
        fxp_rounding = 'floor' if self._half_up else rounding
        self.DATA = Fxp(None, True, dtype=dtype)
        self.DATA.config.rounding = fxp_rounding
        self.DATA.config.overflow = overflow
        # No automatic guard-bit growth: the datapath width is exactly `dtype`
        # (specified via inner_type). DATA_WIDE mirrors DATA so the qw/qcw/shift
        # helpers operate at that same, explicitly-chosen width.
        self.DATA_WIDE = Fxp(None, True, dtype=dtype)
        self.DATA_WIDE.config.rounding = fxp_rounding
        self.DATA_WIDE.config.overflow = overflow
        # half LSB of the target format; exact in float64 for every reachable
        # datapath value (products carry well under 53 significand bits), so
        # the bias itself never adds rounding error
        self._half_lsb = 2.0 ** -(int(self.DATA.n_frac) + 1)

        # fast path: scale to the integer code grid, round, wrap/saturate --
        # exact in float64 (scale is a power of two, codes fit well under 53
        # bits), bit-identical to the fxpmath path above
        fmt = re.fullmatch(r'fxp-s(\d+)/(\d+)', dtype) if isinstance(dtype, str) else None
        self._fast = (fmt is not None and not USE_REFERENCE_FXP
                      and rounding in ('floor', 'around', 'half_up')
                      and overflow in ('wrap', 'saturate'))
        if fmt:
            wbits, fbits = int(fmt.group(1)), int(fmt.group(2))
            self._scale = float(1 << fbits)
            self._inv_scale = 1.0 / self._scale
            self._span = 1 << wbits
            self._max_i = (1 << (wbits - 1)) - 1
            self._min_i = -(1 << (wbits - 1))
        self._wrap = (overflow == 'wrap')
        self._floor = (fxp_rounding == 'floor')

    def _to_int(self, value, half_up_bias):
        x = float(value) * self._scale
        if half_up_bias:
            k = math.floor(x + 0.5)
        elif self._floor:
            k = math.floor(x)
        else:
            k = round(x)            # ties to even, like fxpmath 'around'
        if self._wrap:
            k = ((k - self._min_i) % self._span) + self._min_i
        elif k > self._max_i:
            k = self._max_i
        elif k < self._min_i:
            k = self._min_i
        return k

    def q(self, value):
        if self._fast:
            return self._to_int(value, self._half_up) * self._inv_scale
        if self._half_up:
            value = value + self._half_lsb
        return float(Fxp(value, like=self.DATA))

    def qc(self, value):
        return complex(self.q(np.real(value)), self.q(np.imag(value)))

    def qw(self, value):
        if self._fast:
            return self._to_int(value, self._half_up) * self._inv_scale
        if self._half_up:
            value = value + self._half_lsb
        return float(Fxp(value, like=self.DATA_WIDE))

    def qcw(self, value):
        return complex(self.qw(np.real(value)), self.qw(np.imag(value)))

    def q_shifted(self, value, bits):
        # arithmetic right shift (preadder s0 path and output scaling). In
        # 'half_up' mode the shift rounds: +2^(bits-1) -- half LSB of the
        # result, a carry-in in hardware -- before shifting; otherwise it
        # truncates like the RTL shift_right does today. Note: the pre-shift
        # quantize replicates the reference below, where the bare
        # Fxp(value, dtype=) uses fxpmath DEFAULTS (trunc + saturate), NOT
        # this quantizer's config. Irrelevant in the datapath -- shift inputs
        # are always on-grid and in range -- but kept for bit-identity.
        if bits == 0:
            return self.qw(value)

        if self._fast:
            k = math.trunc(float(value) * self._scale)
            if k > self._max_i:
                k = self._max_i
            elif k < self._min_i:
                k = self._min_i
            if self._half_up:
                k += 1 << (bits - 1)
            k >>= bits              # arithmetic shift, floors like numpy
            return k * self._inv_scale

        fxp_value = Fxp(value, True, dtype=self.DATA_WIDE.dtype)
        shifted = Fxp(None, True, dtype=self.DATA_WIDE.dtype)
        raw = fxp_value.val
        if self._half_up:
            raw = raw + (1 << (bits - 1))
        shifted.val = raw >> bits
        return float(shifted)

    def qcw_shifted(self, value, bits):
        if bits == 0:
            return self.qcw(value)
        return complex(self.q_shifted(np.real(value), bits), self.q_shifted(np.imag(value), bits))


def supported_fft_sizes(limit=3300):
    """The 53 top-level mixed-radix FFT lengths: N = 12 * 2^a * 3^b * 5^c, product <= 275
    (see decomposition_configs.py)."""
    sizes = set()
    for a in range(11):
        for b in range(8):
            for c in range(5):
                p = (2 ** a) * (3 ** b) * (5 ** c)
                if p <= 275 and 12 * p <= limit:
                    sizes.add(12 * p)
    return sorted(sizes)


def smooth_sizes(limit=3300):
    """Every 5-smooth size (2^a * 3^b * 5^c) up to `limit`. Superset of the 53 top-level
    sizes; also covers per-stage sub-sizes (which need not be multiples of 12, e.g. 405).
    Used to fill the increment ROM so any stage/config reconfiguration needs no divide."""
    out = set()
    a = 1
    while a <= limit:
        b = a
        while b <= limit:
            c = b
            while c <= limit:
                out.add(c)
                c *= 5
            b *= 3
        a *= 2
    return sorted(out)


class SinCosLUT:
    """Shared sin/cos table + phase accumulator -- a bit-faithful model of the
    approach-B hardware twiddle generator.

    Hardware modeled (all fixed-point / integer):

      * Phase accumulator, P = addr_bits + frac_bits wide, full circle = 2**P.
        The DDS increment for transform size N is  delta = round(2**P / N)  (one
        divide, done once per reconfiguration). The accumulator value after
        exponent k is  acc = (k * delta) mod 2**P  -- identical to k repeated
        additions of the per-sample increment, because integer add is exact.
          - top `addr_bits` (A)  -> coarse LUT index  (0..L-1),  L = 2**A
          - low  `frac_bits` (f) -> interpolation fraction

      * L/4 quarter-wave ROM: only the first cosine quadrant is stored
        (`cos_q`, L/4+1 words, quantized once with floor rounding like the ROM
        contents). Every other value is produced by EXACT 2's-complement negation
        and address reflection (no re-quantization) -- the standard octant/quadrant
        trick. sin is taken from the same table via sin(t) = cos(t - pi/2).

      * Optional first-order (linear) interpolation between adjacent LUT points,
        with the result requantized to the twiddle word (models the interp
        adder/multiplier output register).

    Error is bounded table + phase quantization; it does NOT accumulate.
    lookup(k, N) returns W_N^k = exp(-j*2*pi*k/N).
    """

    def __init__(self, addr_bits=10, twiddle_dtype='fxp-s16/15', interp=False, frac_bits=14,
                 table_dtype=None, table_rounding='around'):
        self.A = int(addr_bits)
        self.f = int(frac_bits)
        self.L = 1 << self.A                 # circle points
        self.P = self.A + self.f             # phase-accumulator width
        self.Q = self.L >> 2                 # quarter = L/4
        self.interp = bool(interp)
        self.qz = _Quantizer(dtype=twiddle_dtype, rounding=table_rounding)   # interp-output quantizer, twiddles round-to-nearest
        # L/4 ROM trick: store one cosine quadrant only (indices 0..L/4 inclusive).
        # ROM contents are design-time constants, so round-to-nearest (and optionally a
        # wider `table_dtype` than the datapath word) -- this removes the negate-asymmetry
        # and double-quantization loss that floor rounding would add across quadrants.
        qz_tab = _Quantizer(dtype=(table_dtype or twiddle_dtype), rounding=table_rounding)
        self.cos_q = np.array([qz_tab.q(np.cos(2 * np.pi * m / self.L))
                               for m in range(self.Q + 1)])
        self.f_scale = (1.0 / (1 << self.f)) if self.f else 0.0   # constant (interp only)
        self.delta_mem = {}                                       # precomputed: delta_mem[N] = round(2**P/N)

    def precompute(self, sizes):
        """Eagerly fill the universal increment memory `delta_mem[N] = round(2**P/N)`
        for every supported transform (sub-)size. delta depends ONLY on N and P, so this
        one table is identical for every stage -- in hardware it is the same small ROM
        replicated per stage; a stage just addresses the entry for its current sub-size.
        Call once (e.g. lut.precompute(all_sizes)) so no divide ever runs at run time."""
        for N in sizes:
            self.delta_mem[int(N)] = int(round((1 << self.P) / int(N)))
        return self

    def increments(self, N, stage_index, radix):
        """Per-arm increments for a stage = slope * radix**stage_index * delta_mem[N],
        slope = 1..radix-1. `delta_mem[N]` is a memory read (precomputed by precompute());
        the slope / radix**stage_index scalings are small offline constants applied at
        reconfiguration. The run-time datapath then only ADDS these."""
        N = int(N)
        d = self.delta_mem.get(N)
        if d is None:                                    # lazy fallback if not precomputed
            d = int(round((1 << self.P) / N))
            self.delta_mem[N] = d
        base = d * (radix ** stage_index)
        mask = (1 << self.P) - 1
        return [(s * base) & mask for s in range(1, radix)]

    def _sincos(self, p):
        """(cos, sin)(2*pi*p/L) from the single L/4 cosine ROM -- exactly the
        quarter-wave hardware read: index the quarter with `r`, take a second
        (reflected) read `Q-r` for the sine, then apply quadrant sign/swap."""
        p &= (self.L - 1)
        quad = p >> (self.A - 2)          # top 2 address bits -> quadrant
        r = p & (self.Q - 1)              # low A-2 bits -> index into the quarter
        base = self.cos_q[r]              # cos(phi)
        comp = self.cos_q[self.Q - r]     # sin(phi) = cos(pi/2 - phi)  (2nd ROM read)
        if quad == 0:
            return base, comp
        if quad == 1:
            return -comp, base
        if quad == 2:
            return -base, -comp
        return comp, -base

    def lookup_phase(self, acc):
        """P-bit phase-accumulator value -> W = exp(-j*2*pi*acc/2**P) = cos - j*sin."""
        acc &= (1 << self.P) - 1
        if not self.interp:
            i = ((acc + (1 << (self.f - 1))) >> self.f) & (self.L - 1) if self.f else acc & (self.L - 1)
            c, s = self._sincos(i)
            return complex(c, -s)
        i = (acc >> self.f) & (self.L - 1)
        frac = (acc & ((1 << self.f) - 1)) * self.f_scale        # interp fraction
        c0, s0 = self._sincos(i)
        c1, s1 = self._sincos((i + 1) & (self.L - 1))
        cos = self.qz.q(c0 + frac * (c1 - c0))
        sin = self.qz.q(s0 + frac * (s1 - s0))
        return complex(cos, -sin)


class MixedRadix_PreAdder_FXP:
    """Reconfigurable preadder with a CAPABILITY generic (like a VHDL generic /
    if-generate): capability=5 is the full radix-2/3/5 butterfly, capability=3 a
    radix-2/3 one, capability=2 a plain radix-2 butterfly. Only the ports and
    constant multipliers the chosen capability needs are instantiated:

        capability 5: input/output_0..4, constants k2..k6 (4 complex + 1 real mult)
        capability 3: input/output_0..2, constant k6 only   (1 real mult)
        capability 2: input/output_0..1, no multipliers

    The reduced datapaths reuse the exact adder/quantizer sequence of the full
    one's config-2/3 paths, so for the same inputs every capability is bit-exact
    against capability=5."""

    def __init__(self, dtype='fxp-s32/12', shift_bits=0, overflow='saturate', coeff_dtype=None, rounding='floor', capability=5, exit_round=False, internal_frac=None):
        if int(capability) not in [2, 3, 5]:
            raise ValueError('capability must be 5 (radix235), 3 (radix23) or 2 (radix2)')
        self.capability = int(capability)

        self.qz = _Quantizer(dtype=dtype, overflow=overflow, rounding=rounding)
        # Preadder constants k2..k6 live in their own coeff_dtype word. The multiply
        # RESULTS are still quantized through self.qz (the inner_type datapath), so the
        # datapath keeps the data width regardless of the coefficient width.
        self.qzc = _Quantizer(dtype=coeff_dtype, overflow=overflow, rounding=rounding) if coeff_dtype else None
        self.shift_bits = int(shift_bits)
        # exit_round=False (legacy): every internal node re-quantizes at `dtype`
        # (single-dtype convention, bit-exact with today's RTL).
        # exit_round=True: "18 bits in memory, wide in flight" -- inputs are the
        # memory word (`dtype`), the internal adder tree and raw products keep
        # full precision (exact in float64 at these depths; in RTL a wider
        # sfixed on the DSP 25-bit data port / fabric adders), and each output
        # is scaled and rounded ONCE to the memory word before the result
        # FIFOs. `internal_peak` tracks the largest internal magnitude so the
        # required guard bits can be checked against the DSP port budget.
        self.exit_round = bool(exit_round)
        self.internal_peak = 0.0
        # internal_frac (exit_round mode only): fraction bits kept on the raw
        # coefficient products before the t3 adds. None = full precision
        # (data_frac + coeff_frac, the DSP M/P register word); a value like
        # 20-22 models trimming the product so the post-multiplier fabric
        # adders stay narrow. Rounded with the datapath rounding mode. The
        # quantizer's 8 int bits are deliberately generous -- the REAL int
        # requirement comes from internal_peak; this models frac trimming only.
        self.internal_frac = int(internal_frac) if internal_frac else None
        if self.internal_frac and not self.exit_round:
            raise ValueError('internal_frac requires exit_round=True')
        self.qz_prod = (_Quantizer(dtype=f'fxp-s{self.internal_frac + 8}/{self.internal_frac}',
                                   overflow=overflow, rounding=rounding)
                        if self.internal_frac else None)

        self.input_0 = 0.0 + 0.0j
        self.input_1 = 0.0 + 0.0j
        self.output_0 = 0.0 + 0.0j
        self.output_1 = 0.0 + 0.0j
        if self.capability >= 3:
            self.input_2 = 0.0 + 0.0j
            self.output_2 = 0.0 + 0.0j
        if self.capability == 5:
            self.input_3 = 0.0 + 0.0j
            self.input_4 = 0.0 + 0.0j
            self.output_3 = 0.0 + 0.0j
            self.output_4 = 0.0 + 0.0j

        # quantize the constants at coeff_dtype when given; else keep the legacy
        # wide-datapath rounding so existing (single-dtype) runs stay bit-identical.
        qk = self.qzc.q if self.qzc is not None else self.qz.qw
        qkc = self.qzc.qc if self.qzc is not None else self.qz.qcw
        if self.capability == 5:
            self.k2 = qk(0.5 * (np.cos(2 * np.pi / 5) - np.cos(4 * np.pi / 5)))
            self.k3 = qkc(1j * (np.sin(4 * np.pi / 5) - np.sin(2 * np.pi / 5)))
            self.k4 = qkc(-1j * np.sin(4 * np.pi / 5))
            self.k5 = qkc(1j * (np.sin(4 * np.pi / 5) + np.sin(2 * np.pi / 5)))
        if self.capability >= 3:
            self.k6 = qk(-np.sqrt(3) / 2)

    # internal-node quantizers: identity (peak-tracked) in exit_round mode,
    # per-node re-quantization at the datapath word otherwise
    def _qw(self, value):
        if self.exit_round:
            peak = max(abs(value.real), abs(value.imag))
            if peak > self.internal_peak:
                self.internal_peak = peak
            return value
        return self.qz.qcw(value)

    def _qw_shifted(self, value, bits):
        if self.exit_round:
            return self._qw(value * 2.0 ** -bits)   # wide word keeps the bits
        return self.qz.qcw_shifted(value, bits)

    # coefficient-product quantizer: trims the raw product to internal_frac
    # fraction bits (exit_round mode with internal_frac set); otherwise the
    # product follows the internal-node convention (_qw)
    def _qprod(self, value):
        if self.qz_prod is not None:
            return self._qw(self.qz_prod.qcw(value))
        return self._qw(value)

    # output quantizer: scaling shift + ONE rounding to the memory word
    def _qout(self, value):
        if self.shift_bits:
            if self.exit_round:
                return self.qz.qc(value * 2.0 ** -self.shift_bits)
            return self.qz.qc(self.qz.qcw_shifted(value, self.shift_bits))
        return self.qz.qc(value)

    def calculate(self, s0=0, s1=0):
        # capability if-generate: dispatch to the datapath this instance was built with
        if self.capability == 3:
            self._calculate_23(s0)
            return
        if self.capability == 2:
            self._calculate_2()
            return
        # Model the s17/12 input ports of the RTL entity: quantize inputs to the
        # narrow port width before widening (matches resize(i_x*, wide) in VHDL).
        self.input_0 = self.qz.qc(self.input_0)
        self.input_1 = self.qz.qc(self.input_1)
        self.input_2 = self.qz.qc(self.input_2)
        self.input_3 = self.qz.qc(self.input_3)
        self.input_4 = self.qz.qc(self.input_4)

        tmp_0_0 = self._qw(self.input_0)
        tmp_1_0 = self._qw(self.input_1 + self.input_4)
        tmp_2_0 = self._qw(self.input_2 + self.input_3)
        tmp_3_0 = self._qw(self.input_1 - self.input_4)
        tmp_4_0 = self._qw(self.input_2 - self.input_3)

        tmp_0_1 = tmp_0_0
        tmp_1_1 = self._qw(tmp_1_0 + tmp_2_0)
        tmp_2_1 = self._qw(tmp_1_0 - tmp_2_0)
        tmp_3_1 = tmp_3_0
        tmp_4_1 = tmp_4_0
        tmp_5_1 = self._qw(tmp_3_0 + tmp_4_0)

        tmp_0_2 = self._qw(tmp_0_1 + tmp_1_1)

        if s0 == 0:
            tmp_1_2 = self._qw(tmp_0_1 - tmp_1_1)
        elif s0 == 1:
            tmp_1_2 = self._qw(tmp_0_1 - self._qw_shifted(tmp_1_1, 1))
        else:
            tmp_1_2 = self._qw(tmp_0_1 - self._qw_shifted(tmp_1_1, 2))

        mul_1 = self.k6 if (s1 == 0) else self.k2
        mul_1_res = self._qprod(tmp_2_1 * mul_1)

        if s1 == 1:
            tmp_2_2 = mul_1_res
        else:
            tmp_2_2 = self._qw(mul_1_res * 1j)

        tmp_3_2 = self._qprod(tmp_3_1 * self.k3)
        tmp_4_2 = self._qprod(tmp_4_1 * self.k5)
        tmp_5_2 = self._qprod(tmp_5_1 * self.k4)

        tmp_0_3 = tmp_0_2
        tmp_5_3 = tmp_1_2
        tmp_1_3 = self._qw(tmp_1_2 + tmp_2_2)
        tmp_2_3 = self._qw(tmp_1_2 - tmp_2_2)
        tmp_3_3 = self._qw(tmp_3_2 + tmp_5_2)
        tmp_4_3 = self._qw(tmp_4_2 + tmp_5_2)

        # apply optional right-shift scaling before final rounding/quantization
        self.output_0 = self._qout(tmp_0_3)

        tmp_out_1_radix5 = tmp_1_3 + tmp_3_3
        tmp_out_1_radix3 = tmp_1_3
        tmp_out_1_radix2 = tmp_5_3

        if s0 == 0:
            self.output_1 = self._qout(tmp_out_1_radix2)
        elif s0 == 1:
            self.output_1 = self._qout(tmp_out_1_radix3)
        else:
            self.output_1 = self._qout(tmp_out_1_radix5)

        tmp_out_2_radix5 = tmp_2_3 + tmp_4_3
        tmp_out_2_radix3 = tmp_2_3
        out2_val = tmp_out_2_radix3 if (s1 == 0) else tmp_out_2_radix5
        self.output_2 = self._qout(out2_val)

        self.output_4 = self._qout(tmp_1_3 - tmp_3_3)
        self.output_3 = self._qout(tmp_2_3 - tmp_4_3)

    def _calculate_23(self, s0):
        # capability=3 datapath: the full preadder's config-2/3 paths with inputs 3/4
        # tied to zero and the radix-5 multipliers removed. s0: 0 = radix-2, 1 = radix-3.
        self.input_0 = self.qz.qc(self.input_0)
        self.input_1 = self.qz.qc(self.input_1)
        self.input_2 = self.qz.qc(self.input_2)

        tmp_0_0 = self._qw(self.input_0)
        tmp_1_0 = self._qw(self.input_1)
        tmp_2_0 = self._qw(self.input_2)

        tmp_1_1 = self._qw(tmp_1_0 + tmp_2_0)
        tmp_2_1 = self._qw(tmp_1_0 - tmp_2_0)

        tmp_0_2 = self._qw(tmp_0_0 + tmp_1_1)

        if s0 == 0:
            tmp_1_2 = self._qw(tmp_0_0 - tmp_1_1)
        else:
            tmp_1_2 = self._qw(tmp_0_0 - self._qw_shifted(tmp_1_1, 1))

        mul_1_res = self._qprod(tmp_2_1 * self.k6)
        tmp_2_2 = self._qw(mul_1_res * 1j)

        tmp_1_3 = self._qw(tmp_1_2 + tmp_2_2)
        tmp_2_3 = self._qw(tmp_1_2 - tmp_2_2)

        self.output_0 = self._qout(tmp_0_2)

        out1_val = tmp_1_2 if (s0 == 0) else tmp_1_3
        self.output_1 = self._qout(out1_val)

        self.output_2 = self._qout(tmp_2_3)

    def _calculate_2(self):
        # capability=2 datapath: plain 2-point butterfly, no multipliers.
        self.input_0 = self.qz.qc(self.input_0)
        self.input_1 = self.qz.qc(self.input_1)

        tmp_0 = self._qw(self.input_0 + self.input_1)
        tmp_1 = self._qw(self.input_0 - self.input_1)

        self.output_0 = self._qout(tmp_0)
        self.output_1 = self._qout(tmp_1)


class TwiddleROM:
    """Exact reference twiddle source: a precomputed ROM (np.exp), addressed directly by
    the stage's control signals. next_twiddle(phase, delay_cnt) returns
    rom[phase*M + delay_cnt] (arm = phase, position = delay_cnt)."""

    def __init__(self, config, stage_index, size, qz_twiddle):
        self.qzt = qz_twiddle
        self.stage_index = int(stage_index)
        self.reconfigure(config, size)

    def reconfigure(self, config, size):
        R = int(config)
        N = int(size)
        self.M = N // R
        rom = np.ones(N, dtype=complex)
        idx = N - (R - 1) * self.M                # arm-0 block (M samples) stays 1.0
        if self.M > 1:
            for slope in range(1, R):             # arms n1 = 1..R-1
                dk = slope * (R ** self.stage_index)
                k = 0
                for _ in range(self.M):
                    rom[idx] = np.exp(-1j * 2 * np.pi * k / N)
                    k += dk
                    idx += 1
        self.rom = np.array([self.qzt.qc(x) for x in rom], dtype=complex)

    def reset(self):
        pass

    def next_twiddle(self, phase, delay_cnt):
        return self.rom[phase * self.M + delay_cnt]


class TwiddleGenerator:
    """Live twiddle generator (approach B): phase accumulator + shared SinCosLUT, no ROM.

    Driven by the stage's RadixPhaseController signals, passed in each call:
      phase     = arm index n1  (selects which precomputed increment to add)
      delay_cnt = position within the arm-block  (delay_cnt == 0 -> reset accumulator)
    Holds only the accumulator register and the current-increment latch; the per-arm
    increments are precomputed constants (LUT increment memory). The datapath only ADDS."""

    def __init__(self, config, stage_index, size, lut, qz_twiddle, configs=None):
        self.lut = lut
        self.pmask = (1 << lut.P) - 1
        self.one = qz_twiddle.qc(1 + 0j)
        self.stage_index = int(stage_index)         # fixed per physical stage
        # ---- FIXED increment ROM ----
        # Precompute, ONCE, the per-arm phase increments for every (radix, size) this
        # stage can be configured to. All divides / multiplies / powers live HERE (offline
        # ROM generation). Reconfiguring the radix is then a pure lookup -- no arithmetic.
        # Each entry: (arm_inc list, rotate_on flag).
        if configs is None:
            configs = [(R, N) for N in smooth_sizes() for R in (2, 3, 5) if N % R == 0]
        lut.precompute(sorted({N for (_, N) in configs}))          # delta_mem: the divides
        # Flat increment ROM: each config's (R-1) arm increments occupy a contiguous
        # block. A per-config BASE ADDRESS points at the block; the arm index (n1)
        # offsets into it, so one increment is read at a time (per arm), not all at once.
        self.inc_rom = []                          # linear ROM of increment words
        self._base = {}                            # config (R,N) -> base address
        self._rot_on = {}                          # config (R,N) -> rotation active?
        for (R, N) in configs:
            self._base[(R, N)] = len(self.inc_rom)
            d = self.lut.delta_mem[N] * (R ** self.stage_index)     # offline
            for s in range(1, R):
                self.inc_rom.append((s * d) & self.pmask)           # offline
            self._rot_on[(R, N)] = (N // R) > 1
        self.reconfigure(config, size)

    def reconfigure(self, config, size):
        """Switch radix/size: look up just the ROM BASE ADDRESS for this config (an
        offset, not the data). Increments are read one-by-one from inc_rom during run."""
        self.R = int(config)
        key = (self.R, int(size))
        self.base = self._base[key]                # base address (offset) into inc_rom
        self.rotate_on = self._rot_on[key]
        self.reset()

    def reset(self):
        self.acc = 0
        self.cur_inc = 0

    def next_twiddle(self, phase, delay_cnt):
        if delay_cnt == 0:                         # arm boundary: read THIS arm's increment
            self.acc = 0
            self.cur_inc = self.inc_rom[self.base + phase - 1] if phase >= 1 else 0
        if (not self.rotate_on) or phase == 0:
            return self.one                        # arm-0 region: no rotation (W = 1)
        tw = self.lut.lookup_phase(self.acc)               # sin/cos ROM read
        self.acc = (self.acc + self.cur_inc) & self.pmask  # accumulator ADD
        return tw


class TwiddleOctantROM:
    """Stored-twiddle source with octant/quarter/half symmetry, holding ALL of a stage's
    configs in ONE flat ROM. `reconfigure(config, size)` SELECTS that config's block
    (base address / N / reduction level) -- it does NOT rebuild, so at run time you can
    switch (radix, N) and the twiddles are already resident (hardware-faithful).

    Per config the reduced set W_N^n is stored (quantized to twiddle_dtype); any W_N^k is
    rebuilt by exact conj/negate/swap (no multiplier). Reduction level from N:
        N % 4 == 0 -> octant  (N/8 + 1) ;  N even -> quarter (N/4 + 1) ;  N odd -> half."""

    def __init__(self, configs, stage_index, qz_twiddle):
        self.stage_index = int(stage_index)
        self.one = qz_twiddle.qc(1 + 0j)
        self.rom = []                 # flat ROM: every config's reduced block, concatenated
        self.meta = {}                # (R, N) -> dict(base, N, level, R, rotate_on)
        for (R, N) in configs:
            R = int(R); N = int(N)
            if (R, N) in self.meta:
                continue
            if N % 4 == 0:
                lvl, bound = 3, N // 8       # octant
            elif N % 2 == 0:
                lvl, bound = 2, N // 4       # quarter
            else:
                lvl, bound = 1, (N - 1) // 2  # half (conjugate only)
            base = len(self.rom)
            for n in range(bound + 1):
                self.rom.append(qz_twiddle.qc(np.exp(-1j * 2 * np.pi * n / N)))
            self.meta[(R, N)] = dict(base=base, N=N, level=lvl, R=R, rotate_on=(N // R) > 1)
        self.cur = None
        if configs:
            self.reconfigure(*configs[0])

    def reconfigure(self, config, size):
        self.cur = self.meta[(int(config), int(size))]   # pure select -- no rebuild

    def reset(self):
        pass

    def _recon(self, k, c):
        """Reconstruct W_N^k from config `c`'s reduced block by exact symmetry folds.
        Transforms are wiring + negate in hardware (no multiplier):
          conj  = negate Im ;  x(-1) = negate both ;  x(-j) = swap Re/Im, negate new Im."""
        N = c['N']; base = c['base']; level = c['level']
        k %= N
        m = 1 + 0j
        cj = False

        def apply(t, cc):
            nonlocal m, cj
            m = m * (np.conj(t) if cj else t)
            cj = cj != cc

        if 2 * k > N:                    # conjugate:  W_k = conj(W_{N-k})
            k = N - k;        apply(1, True)
        if level >= 2 and 4 * k > N:     # quarter:    W_k = -conj(W_{N/2-k})
            k = (N // 2) - k; apply(-1, True)
        if level >= 3 and 8 * k > N:     # octant:     W_k = -j*conj(W_{N/4-k})
            k = (N // 4) - k; apply(-1j, True)
        v = self.rom[base + k]           # base(config) + reduced index
        if cj:
            v = complex(v.real, -v.imag)
        return complex(m * v)

    def next_twiddle(self, phase, delay_cnt):
        c = self.cur
        if (not c['rotate_on']) or phase == 0:
            return self.one
        k = phase * delay_cnt * (c['R'] ** self.stage_index)   # twiddle exponent
        return self._recon(k, c)


class MixedRadix_Rotator_FXP:
    """Rotator: output = input * twiddle. The twiddle is produced by a separate source
    component -- TwiddleGenerator (LUT, 'lut'), TwiddleOctantROM (stored, 'octant'), or
    TwiddleROM (exact reference, default). The rotator itself only rotates (one mult)."""

    def __init__(self, config, stage_index, size, dtype='fxp-s32/12', overflow='saturate', twiddle_dtype=None,
                 twiddle_source='exact', lut=None, twiddle_rounding='around', twiddle_configs=None, rounding='floor'):
        if int(config) not in [2, 3, 5]:
            raise ValueError('config must be 2, 3 or 5')

        self.qz = _Quantizer(dtype=dtype, overflow=overflow, rounding=rounding)   # data multiply quantizer
        if twiddle_dtype is None:
            twiddle_dtype = dtype
        # Twiddle values may round-to-nearest, independent of the floor-rounded data path.
        qz_twiddle = _Quantizer(dtype=twiddle_dtype, overflow=overflow, rounding=twiddle_rounding)

        # twiddle_configs: the full (radix, size) set this stage must support. If given,
        # the ROM holds all of them and reconfigure() just selects; else it defaults to
        # the single (config, size) built here.
        cfgs = [(int(r), int(n)) for (r, n) in twiddle_configs] if twiddle_configs else [(int(config), int(size))]

        src = str(twiddle_source)
        if src == 'lut' and lut is not None:
            self.twiddle = TwiddleGenerator(config, stage_index, size, lut, qz_twiddle,
                                            configs=(cfgs if twiddle_configs else None))
        elif src == 'octant':
            self.twiddle = TwiddleOctantROM(cfgs, stage_index, qz_twiddle)
        else:
            self.twiddle = TwiddleROM(config, stage_index, size, qz_twiddle)

        self.config = int(config)
        self.size = int(size)
        self.input = 0.0 + 0.0j
        self.output = 0.0 + 0.0j

    def reset(self):
        self.twiddle.reset()

    def reconfigure(self, config, size):
        self.twiddle.reconfigure(config, size)
        self.config = int(config)
        self.size = int(size)

    def rotate(self, fifo_full_flag, phase, delay_cnt):
        if fifo_full_flag:
            self.output = self.qz.qc(self.input * self.twiddle.next_twiddle(phase, delay_cnt))


class RadixPhaseController:
    def __init__(self, cfg, cfg_delay):
        self.cfg = int(cfg)
        self.cfg_delay = int(cfg_delay)
        if self.cfg_delay < 1:
            raise ValueError('cfg_delay must be >= 1')
        self.phase = 0
        self.delay_cnt = 0

    def reset(self):
        self.phase = 0
        self.delay_cnt = 0

    def tick(self, valid=True):
        if not valid:
            return
        if self.delay_cnt == (self.cfg_delay - 1):
            self.delay_cnt = 0
            self.phase = 0 if self.phase == (self.cfg - 1) else self.phase + 1
        else:
            self.delay_cnt += 1


class MixedRadix_SDF_stage_counter_ctrl_FXP:
    """SDF butterfly stage with a CAPABILITY generic (mirrors the stage types of
    docs/reconfigurable_fft_midbypass_reversed_architecture.md):

        capability 5 (radix235): 4 delay FIFOs, full preadder, config in {2,3,5}
        capability 3 (radix23) : 2 delay FIFOs, k6-only preadder, config in {2,3}
        capability 2 (radix2)  : 1 delay FIFO, adder-only preadder, config = 2

    The runtime `config` (radix) must be within the build-time capability. For the
    same (config, size) every capability produces a bit-identical output stream, so
    stages of different capability can be mixed freely in a chain."""

    def __init__(self, config, stage_index, size, cfg_delay=None, s0=2, s1=1, dtype='fxp-s32/12', preadder_shift_bits=None, overflow='wrap', twiddle_dtype=None, twiddle_source='exact', lut=None, twiddle_rounding='around', twiddle_configs=None, inner_type=None, input_dtype=None, coeff_dtype=None, rounding='floor', capability=5, preadder_exit_round=False, preadder_internal_frac=None):
        if int(capability) not in [2, 3, 5]:
            raise ValueError('capability must be 5 (radix235), 3 (radix23) or 2 (radix2)')
        self.capability = int(capability)
        if int(config) not in [2, 3, 5]:
            raise ValueError('This implementation is for radix-2-3-5 stages (config=5 or 3 or 2)')
        if int(config) > self.capability:
            raise ValueError(f'config={int(config)} exceeds stage capability {self.capability}')

        # `inner_type` is the datapath word for this stage (falls back to `dtype`).
        inner = inner_type or dtype
        self.qz = _Quantizer(dtype=inner, overflow=overflow, rounding=rounding)
        # `input_dtype`, when set (only the first stage of a chain does), quantizes the
        # incoming sample to the chain input word before it enters the inner datapath.
        self.qz_in = _Quantizer(dtype=input_dtype, overflow=overflow, rounding=rounding) if input_dtype else None
        self.config = int(config)
        self.stage_index = int(stage_index)
        self.num_of_samples = int(size)

        if cfg_delay is None:
            if self.config == 5:
                self.cfg_delay = self.num_of_samples // 5
            elif self.config == 3:
                self.cfg_delay = self.num_of_samples // 3
            else:
                self.cfg_delay = self.num_of_samples // 2
        else:
            self.cfg_delay = int(cfg_delay)

        if self.cfg_delay < 1:
            raise ValueError('cfg_delay must be >= 1')

        self.s0 = int(s0)
        self.s1 = int(s1)

        # capability if-generate: instantiate only the FIFOs this stage type has
        self.fifo_0 = Fifo(self.cfg_delay)
        if self.capability >= 3:
            self.fifo_1 = Fifo(self.cfg_delay)
        if self.capability == 5:
            self.fifo_2 = Fifo(self.cfg_delay)
            self.fifo_3 = Fifo(self.cfg_delay)

        if preadder_shift_bits is None:
            if self.config == 2:
                preadder_shift_bits = 1
            elif self.config == 3:
                preadder_shift_bits = 2
            else:
                preadder_shift_bits = 3

        self.pre_adder = MixedRadix_PreAdder_FXP(dtype=inner, shift_bits=preadder_shift_bits, overflow=overflow, coeff_dtype=coeff_dtype, rounding=rounding, capability=self.capability, exit_round=preadder_exit_round, internal_frac=preadder_internal_frac)
        self.rotator = MixedRadix_Rotator_FXP(config=self.config, stage_index=self.stage_index, size=self.num_of_samples, dtype=inner, overflow=overflow, twiddle_dtype=twiddle_dtype, twiddle_source=twiddle_source, lut=lut, twiddle_rounding=twiddle_rounding, twiddle_configs=twiddle_configs, rounding=rounding)

        self.ctrl = RadixPhaseController(cfg=self.config, cfg_delay=self.cfg_delay)
        # Twiddle control is just self.ctrl DELAYED to the twiddle-multiply pipeline
        # stage: push (phase, delay_cnt) every cycle, pop when the rotator rotates
        # (rot_en). This delay line auto-matches the stage's datapath latency, so the
        # generator always sees the control that belongs to the data being rotated.
        self._tw_delay = deque()

        self.input_sample = 0.0 + 0.0j
        self.output_sample = 0.0 + 0.0j
        self.op_cnt = 0

    def set_muxes(self, s0, s1):
        self.s0 = int(s0)
        self.s1 = int(s1)

    def reset(self):
        self.ctrl.reset()
        self._tw_delay.clear()
        self.rotator.reset()
        self.input_sample = 0.0 + 0.0j
        self.output_sample = 0.0 + 0.0j
        self.op_cnt = 0

    def _load_preadder_inputs(self):
        # drive only the ports the preadder capability instantiates; unused ports of a
        # WIDER capability are tied to zero (exactly the original radix235 behavior)
        self.pre_adder.input_0 = self.qz.qc(self.fifo_0.get_output())

        if self.config == 5:
            self.pre_adder.input_1 = self.qz.qc(self.fifo_1.get_output())
            self.pre_adder.input_2 = self.qz.qc(self.fifo_2.get_output())
            self.pre_adder.input_3 = self.qz.qc(self.fifo_3.get_output())
            self.pre_adder.input_4 = self.input_sample
        elif self.config == 3:
            self.pre_adder.input_1 = self.qz.qc(self.fifo_1.get_output())
            self.pre_adder.input_2 = self.input_sample
            if self.capability == 5:
                self.pre_adder.input_3 = 0.0 + 0.0j
                self.pre_adder.input_4 = 0.0 + 0.0j
        else:
            self.pre_adder.input_1 = self.input_sample
            if self.capability >= 3:
                self.pre_adder.input_2 = 0.0 + 0.0j
            if self.capability == 5:
                self.pre_adder.input_3 = 0.0 + 0.0j
                self.pre_adder.input_4 = 0.0 + 0.0j

    def calculate(self, input_sample, valid=True):
        if not valid:
            return self.output_sample

        # chain input port (first stage) -> input_dtype, then the inner datapath word
        if self.qz_in is not None:
            self.input_sample = self.qz.qc(self.qz_in.qc(input_sample))
        else:
            self.input_sample = self.qz.qc(input_sample)

        self._load_preadder_inputs()

        if self.config == 2:
            self.s0 = 0
        elif self.config == 3:
            self.s0 = 1
            self.s1 = 0
        else:
            self.s0 = 2
            self.s1 = 1

        self.pre_adder.calculate(s0=self.s0, s1=self.s1)

        phase = self.ctrl.phase

        if phase == 0:
            pre_rot = self.qz.qc(self.fifo_0.get_output())
            self.fifo_0.shift(self.input_sample)
        elif phase == 1:
            pre_rot = self.qz.qc(self.fifo_1.get_output()) if (self.config != 2) else self.pre_adder.output_0
            if self.config != 2:
                self.fifo_1.shift(self.input_sample)
            else:
                self.fifo_0.shift(self.pre_adder.output_1)
        elif phase == 2:
            pre_rot = self.qz.qc(self.fifo_2.get_output()) if (self.config != 3) else self.pre_adder.output_0
            if self.config != 3:
                self.fifo_2.shift(self.input_sample)
            else:
                self.fifo_0.shift(self.pre_adder.output_1)
                self.fifo_1.shift(self.pre_adder.output_2)
        elif phase == 3:
            pre_rot = self.qz.qc(self.fifo_3.get_output())
            self.fifo_3.shift(self.input_sample)
        else:
            pre_rot = self.pre_adder.output_0
            self.fifo_0.shift(self.pre_adder.output_1)
            self.fifo_1.shift(self.pre_adder.output_2)
            self.fifo_2.shift(self.pre_adder.output_3)
            self.fifo_3.shift(self.pre_adder.output_4)

        self.rotator.input = self.qz.qc(pre_rot)

        if self.config == 2:
            rot_en = self.fifo_0.is_full()
        elif self.config == 3:
            rot_en = self.fifo_1.is_full()
        else:
            rot_en = self.fifo_3.is_full()

        # delay the phase controller signals to the twiddle-multiply stage
        self._tw_delay.append((self.ctrl.phase, self.ctrl.delay_cnt))
        if rot_en:
            tw_phase, tw_dc = self._tw_delay.popleft()
        else:
            tw_phase, tw_dc = 0, 0        # not used (rotate is a no-op when not rot_en)
        self.rotator.rotate(rot_en, tw_phase, tw_dc)
        self.output_sample = self.qz.qc(self.rotator.output)

        self.ctrl.tick(valid=True)

        if self.op_cnt == (self.num_of_samples - 1):
            self.op_cnt = 0
        else:
            self.op_cnt += 1

        return self.output_sample


class MixedRadix_FinalScaler_FXP:
    def __init__(self, size, total_shift=0, dtype='fxp-s32/12', overflow='saturate', inner_type=None, output_dtype=None, rounding='floor'):
        self.size = int(size)
        self.total_shift = int(total_shift)
        # input arrives on the inner datapath word; the scaled result is the chain
        # output word (output_dtype). Both fall back to `dtype`.
        inner = inner_type or dtype
        out = output_dtype or dtype
        self.qz = _Quantizer(dtype=inner, overflow=overflow, rounding=rounding)
        self.qz_out = _Quantizer(dtype=out, overflow=overflow, rounding=rounding)
        self.input = 0.0 + 0.0j
        self.output = 0.0 + 0.0j

        if self.size < 1:
            raise ValueError('size must be >= 1')

        self.scale = (2 ** self.total_shift) / self.size

    def scale_sample(self, value):
        self.input = self.qz.qc(value)
        self.output = self.qz_out.qc(self.input * self.scale)
        # print("SCALE = ", self.scale)
        return self.output
