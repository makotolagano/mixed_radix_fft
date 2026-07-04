import re

import numpy as np
from fxpmath import Fxp
from utils import Fifo


class _Quantizer:
    def __init__(self, dtype='fxp-s32/12', internal_growth_bits=3, overflow='wrap'):
        self.DATA = Fxp(None, True, dtype=dtype)
        self.DATA.config.rounding = 'floor'
        self.DATA.config.overflow = overflow
        self.DATA_WIDE = Fxp(None, True, dtype=self._grow_dtype(dtype, internal_growth_bits))
        self.DATA_WIDE.config.rounding = 'floor'
        self.DATA_WIDE.config.overflow = overflow

    @staticmethod
    def _grow_dtype(dtype, internal_growth_bits):
        match = re.fullmatch(r'fxp-s(\d+)/(\d+)', str(dtype))
        if match is None:
            return dtype

        total_bits = int(match.group(1)) + int(internal_growth_bits)
        frac_bits = int(match.group(2)) + int(internal_growth_bits)
        return f'fxp-s{total_bits}/{frac_bits}'

    def q(self, value):
        return float(Fxp(value, like=self.DATA))

    def qc(self, value):
        return complex(self.q(np.real(value)), self.q(np.imag(value)))

    def qw(self, value):
        return float(Fxp(value, like=self.DATA_WIDE))

    def qcw(self, value):
        return complex(self.qw(np.real(value)), self.qw(np.imag(value)))

    def q_shifted(self, value, bits):
        if bits == 0:
            return self.qw(value)

        fxp_value = Fxp(value, True, dtype=self.DATA_WIDE.dtype)
        shifted = Fxp(None, True, dtype=self.DATA_WIDE.dtype)
        shifted.val = fxp_value.val >> bits
        return float(shifted)

    def qcw_shifted(self, value, bits):
        if bits == 0:
            return self.qcw(value)
        return complex(self.q_shifted(np.real(value), bits), self.q_shifted(np.imag(value), bits))


class MixedRadix_PreAdder_FXP:
    def __init__(self, dtype='fxp-s32/12', shift_bits=0, overflow='wrap'):
        self.qz = _Quantizer(dtype=dtype, overflow=overflow)
        self.shift_bits = int(shift_bits)

        self.input_0 = 0.0 + 0.0j
        self.input_1 = 0.0 + 0.0j
        self.input_2 = 0.0 + 0.0j
        self.input_3 = 0.0 + 0.0j
        self.input_4 = 0.0 + 0.0j

        self.output_0 = 0.0 + 0.0j
        self.output_1 = 0.0 + 0.0j
        self.output_2 = 0.0 + 0.0j
        self.output_3 = 0.0 + 0.0j
        self.output_4 = 0.0 + 0.0j

        self.k2 = self.qz.qw(0.5 * (np.cos(2 * np.pi / 5) - np.cos(4 * np.pi / 5)))
        self.k3 = self.qz.qcw(1j * (np.sin(4 * np.pi / 5) - np.sin(2 * np.pi / 5)))
        self.k4 = self.qz.qcw(-1j * np.sin(4 * np.pi / 5))
        self.k5 = self.qz.qcw(1j * (np.sin(4 * np.pi / 5) + np.sin(2 * np.pi / 5)))
        self.k6 = self.qz.qw(-np.sqrt(3) / 2)

    def calculate(self, s0, s1):
        self.input_0 = self.qz.qc(self.input_0)
        self.input_1 = self.qz.qc(self.input_1)
        self.input_2 = self.qz.qc(self.input_2)
        self.input_3 = self.qz.qc(self.input_3)
        self.input_4 = self.qz.qc(self.input_4)

        tmp_0_0 = self.qz.qcw(self.input_0)
        tmp_1_0 = self.qz.qcw(self.input_1 + self.input_4)
        tmp_2_0 = self.qz.qcw(self.input_2 + self.input_3)
        tmp_3_0 = self.qz.qcw(self.input_1 - self.input_4)
        tmp_4_0 = self.qz.qcw(self.input_2 - self.input_3)

        tmp_0_1 = tmp_0_0
        tmp_1_1 = self.qz.qcw(tmp_1_0 + tmp_2_0)
        tmp_2_1 = self.qz.qcw(tmp_1_0 - tmp_2_0)
        tmp_3_1 = tmp_3_0
        tmp_4_1 = tmp_4_0
        tmp_5_1 = self.qz.qcw(tmp_3_0 + tmp_4_0)

        tmp_0_2 = self.qz.qcw(tmp_0_1 + tmp_1_1)

        if s0 == 0:
            tmp_1_2 = self.qz.qcw(tmp_0_1 - tmp_1_1)
        elif s0 == 1:
            tmp_1_2 = self.qz.qcw(tmp_0_1 - self.qz.qcw_shifted(tmp_1_1, 1))
        else:
            tmp_1_2 = self.qz.qcw(tmp_0_1 - self.qz.qcw_shifted(tmp_1_1, 2))

        mul_1 = self.k6 if (s1 == 0) else self.k2
        mul_1_res = self.qz.qcw(tmp_2_1 * mul_1)

        if s1 == 1:
            tmp_2_2 = mul_1_res
        else:
            tmp_2_2 = self.qz.qcw(mul_1_res * 1j)

        tmp_3_2 = self.qz.qcw(tmp_3_1 * self.k3)
        tmp_4_2 = self.qz.qcw(tmp_4_1 * self.k5)
        tmp_5_2 = self.qz.qcw(tmp_5_1 * self.k4)

        tmp_0_3 = tmp_0_2
        tmp_5_3 = tmp_1_2
        tmp_1_3 = self.qz.qcw(tmp_1_2 + tmp_2_2)
        tmp_2_3 = self.qz.qcw(tmp_1_2 - tmp_2_2)
        tmp_3_3 = self.qz.qcw(tmp_3_2 + tmp_5_2)
        tmp_4_3 = self.qz.qcw(tmp_4_2 + tmp_5_2)

        # apply optional right-shift scaling before final rounding/quantization
        self.output_0 = self.qz.qc(self.qz.qcw_shifted(tmp_0_3, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_0_3)

        tmp_out_1_radix5 = tmp_1_3 + tmp_3_3
        tmp_out_1_radix3 = tmp_1_3
        tmp_out_1_radix2 = tmp_5_3

        if s0 == 0:
            self.output_1 = self.qz.qc(self.qz.qcw_shifted(tmp_out_1_radix2, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_out_1_radix2)
        elif s0 == 1:
            self.output_1 = self.qz.qc(self.qz.qcw_shifted(tmp_out_1_radix3, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_out_1_radix3)
        else:
            self.output_1 = self.qz.qc(self.qz.qcw_shifted(tmp_out_1_radix5, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_out_1_radix5)

        tmp_out_2_radix5 = tmp_2_3 + tmp_4_3
        tmp_out_2_radix3 = tmp_2_3
        out2_val = tmp_out_2_radix3 if (s1 == 0) else tmp_out_2_radix5
        self.output_2 = self.qz.qc(self.qz.qcw_shifted(out2_val, self.shift_bits)) if self.shift_bits else self.qz.qc(out2_val)

        self.output_4 = self.qz.qc(self.qz.qcw_shifted(tmp_1_3 - tmp_3_3, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_1_3 - tmp_3_3)
        self.output_3 = self.qz.qc(self.qz.qcw_shifted(tmp_2_3 - tmp_4_3, self.shift_bits)) if self.shift_bits else self.qz.qc(tmp_2_3 - tmp_4_3)


class MixedRadix_Rotator_FXP:
    cnt = 0

    def __init__(self, config, stage_index, size, dtype='fxp-s32/12', overflow='saturate', twiddle_dtype=None):
        if int(config) not in [2, 3, 5]:
            raise ValueError('config must be 2, 3 or 5')

        self.qz = _Quantizer(dtype=dtype, overflow=overflow)
        # Separate quantizer for twiddle factors; defaults to same as data if not specified
        if twiddle_dtype is None:
            twiddle_dtype = dtype
        self.qz_twiddle = _Quantizer(dtype=twiddle_dtype, overflow=overflow)
        self.config = int(config)
        self.stage_index = int(stage_index)
        self.size = int(size)
        self.input = 0.0 + 0.0j
        self.output = 0.0 + 0.0j

        self.twiddleROM = np.ones(self.size, dtype=complex)
        N = self.size

        if self.config == 2:
            active_len = self.size // 2
            if active_len > 1:
                for i in range(active_len):
                    k = i * 2 ** self.stage_index
                    self.twiddleROM[i + (self.size - active_len)] = np.exp(-1j * 2 * np.pi * k / N)
        elif self.config == 3:
            active_len = self.size - self.size // 3
            if active_len > 2:
                for i in range(active_len):
                    if i < active_len // 2:
                        k = i * 3 ** self.stage_index
                    else:
                        k = 2 * (i - active_len // 2) * 3 ** self.stage_index
                    self.twiddleROM[i + (self.size - active_len)] = np.exp(-1j * 2 * np.pi * k / N)
        else:
            active_len = self.size - self.size // 5
            if active_len > 4:
                for i in range(active_len):
                    if i < active_len / 4:
                        k = i * 5 ** self.stage_index
                    elif i < active_len / 2:
                        k = 2 * (i - active_len // 4) * 5 ** self.stage_index
                    elif i < 3 * active_len / 4:
                        k = 3 * (i - 2 * active_len // 4) * 5 ** self.stage_index
                    else:
                        k = 4 * (i - 3 * active_len // 4) * 5 ** self.stage_index
                    self.twiddleROM[i + (self.size - active_len)] = np.exp(-1j * 2 * np.pi * k / N)

        self.twiddleROM = np.array([self.qz_twiddle.qc(x) for x in self.twiddleROM], dtype=complex)

    def rotate(self, fifo_full_flag):
        if fifo_full_flag:
            self.output = self.qz.qc(self.input * self.twiddleROM[self.cnt])
            if self.cnt == len(self.twiddleROM) - 1:
                self.cnt = 0
            else:
                self.cnt += 1


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
    def __init__(self, config, stage_index, size, cfg_delay=None, s0=2, s1=1, dtype='fxp-s32/12', preadder_shift_bits=None, overflow='saturate', twiddle_dtype=None):
        if int(config) not in [2, 3, 5]:
            raise ValueError('This implementation is for radix-2-3-5 stages (config=5 or 3 or 2)')

        self.qz = _Quantizer(dtype=dtype, overflow=overflow)
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

        self.fifo_0 = Fifo(self.cfg_delay)
        self.fifo_1 = Fifo(self.cfg_delay)
        self.fifo_2 = Fifo(self.cfg_delay)
        self.fifo_3 = Fifo(self.cfg_delay)

        if preadder_shift_bits is None:
            if self.config == 2:
                preadder_shift_bits = 1
            elif self.config == 3:
                preadder_shift_bits = 2
            else:
                preadder_shift_bits = 3

        self.pre_adder = MixedRadix_PreAdder_FXP(dtype=dtype, shift_bits=preadder_shift_bits, overflow=overflow)
        self.rotator = MixedRadix_Rotator_FXP(config=self.config, stage_index=self.stage_index, size=self.num_of_samples, dtype=dtype, overflow=overflow, twiddle_dtype=twiddle_dtype)

        self.ctrl = RadixPhaseController(cfg=self.config, cfg_delay=self.cfg_delay)

        self.input_sample = 0.0 + 0.0j
        self.output_sample = 0.0 + 0.0j
        self.op_cnt = 0

    def set_muxes(self, s0, s1):
        self.s0 = int(s0)
        self.s1 = int(s1)

    def reset(self):
        self.ctrl.reset()
        self.input_sample = 0.0 + 0.0j
        self.output_sample = 0.0 + 0.0j
        self.op_cnt = 0

    def _load_preadder_inputs(self):
        self.pre_adder.input_0 = self.qz.qc(self.fifo_0.get_output())
        self.pre_adder.input_1 = self.qz.qc(self.fifo_1.get_output()) if self.config != 2 else self.input_sample

        if self.config == 5:
            self.pre_adder.input_2 = self.qz.qc(self.fifo_2.get_output())
        elif self.config == 3:
            self.pre_adder.input_2 = self.input_sample
        else:
            self.pre_adder.input_2 = 0.0 + 0.0j

        self.pre_adder.input_3 = self.qz.qc(self.fifo_3.get_output()) if self.config == 5 else 0.0 + 0.0j
        self.pre_adder.input_4 = self.input_sample if self.config == 5 else 0.0 + 0.0j

    def calculate(self, input_sample, valid=True):
        if not valid:
            return self.output_sample

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

        self.rotator.rotate(rot_en)
        self.output_sample = self.qz.qc(self.rotator.output)

        self.ctrl.tick(valid=True)

        if self.op_cnt == (self.num_of_samples - 1):
            self.op_cnt = 0
        else:
            self.op_cnt += 1

        return self.output_sample


class MixedRadix_FinalScaler_FXP:
    def __init__(self, size, total_shift=0, dtype='fxp-s32/12', internal_growth_bits=3, overflow='saturate'):
        self.size = int(size)
        self.total_shift = int(total_shift)
        self.qz = _Quantizer(dtype=dtype, internal_growth_bits=internal_growth_bits, overflow=overflow)
        self.input = 0.0 + 0.0j
        self.output = 0.0 + 0.0j

        if self.size < 1:
            raise ValueError('size must be >= 1')

        self.scale = (2 ** self.total_shift) / self.size

    def scale_sample(self, value):
        self.input = self.qz.qc(value)
        self.output = self.qz.qc(self.input * self.scale)
        # print("SCALE = ", self.scale)
        return self.output
