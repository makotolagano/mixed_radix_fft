import numpy as np
from scipy.fft import fft
from utils import Fifo

class MixedRadix_PreAdder:
    def __init__(self):
        self.input_0 = 0.0
        self.input_1 = 0.0
        self.input_2 = 0.0
        self.input_3 = 0.0
        self.input_4 = 0.0
        self.output_0 = 0.0
        self.output_1 = 0.0
        self.output_2 = 0.0
        self.output_3 = 0.0
        self.output_4 = 0.0

        self.k1 = -1/4
        self.k2 = 1/2 * (np.cos(2*np.pi/5) - np.cos(4*np.pi/5))
        self.k3 = 1j * (np.sin(4*np.pi/5) - np.sin(2*np.pi/5))
        self.k4 = -1j * np.sin(4*np.pi/5)
        self.k5 = 1j * (np.sin(4*np.pi/5) + np.sin(2*np.pi/5))

        self.k6 = -np.sqrt(3)/2


    def calculate(self, s0, s1):
        tmp_0_0 = self.input_0
        tmp_1_0 = self.input_1 + self.input_4
        tmp_2_0 = self.input_2 + self.input_3
        tmp_3_0 = self.input_1 - self.input_4
        tmp_4_0 = self.input_2 - self.input_3
        ###### prvi nivo pajplajna ^
        tmp_0_1 = tmp_0_0
        tmp_1_1 = tmp_1_0 + tmp_2_0
        tmp_2_1 = tmp_1_0 - tmp_2_0
        tmp_3_1 = tmp_3_0
        tmp_4_1 = tmp_4_0
        tmp_5_1 = tmp_3_0 + tmp_4_0 # dodatna grana izmedju
        ###### drugi nivo pajplajna ^
        tmp_0_2 = tmp_0_1 + tmp_1_1

        mul_0 = 0.0
        if (s0 == 0):
            mul_0 = -1.0
        elif (s0 == 1):
            mul_0 = -0.5
        elif (s0 == 2):
            mul_0 = -0.25

        tmp_1_2 = tmp_0_1 + tmp_1_1 * mul_0 #(-0.25 or -0.5 or -1.0)
        # print("tmp_1_2 = ", tmp_1_2.real, " +j ", tmp_1_2.imag)

        mul_1 = 0.0
        if (s1 == 0):
            mul_1 = self.k6
        elif (s1 == 1):
            mul_1 = self.k2
        
        mul_1_res = tmp_2_1 * mul_1

        tmp_2_2 = mul_1_res if (s1 == 1) else mul_1_res * (1j)
        # print("tmp_2_2 = ", tmp_2_2.real, " +j ", tmp_2_2.imag)

        tmp_3_2 = tmp_3_1 * self.k3 #(-1j*0.363)
        # print("tmp_3_2 = ", tmp_3_2.real, " +j ", tmp_3_2.imag)
        tmp_4_2 = tmp_4_1 * self.k5 #1j*1.539
        # print("tmp_4_2 = ", tmp_4_2.real, " +j ", tmp_4_2.imag)
        tmp_5_2 = tmp_5_1 * self.k4 #(-1j*0.588)
        # print("tmp_5_2 = ", tmp_5_2.real, " +j ", tmp_5_2.imag)
        ###### treci nivo pajplajna ^
        tmp_0_3 = tmp_0_2
        tmp_5_3 = tmp_1_2
        tmp_1_3 = tmp_1_2 + tmp_2_2
        tmp_2_3 = tmp_1_2 - tmp_2_2
        tmp_3_3 = tmp_3_2 + tmp_5_2
        tmp_4_3 = tmp_4_2 + tmp_5_2
        ###### cetvrti nivo pajplajna ^
        self.output_0 = tmp_0_3

        tmp_out_1_radix5 = tmp_1_3 + tmp_3_3
        tmp_out_1_radix3 = tmp_1_3
        tmp_out_1_radix2 = tmp_5_3
        
        if (s0 == 0):
            self.output_1 = tmp_out_1_radix2
        elif (s0 == 1):
            self.output_1 = tmp_out_1_radix3
        elif (s0 == 2):
            self.output_1 = tmp_out_1_radix5
        
        tmp_out_2_radix5 = tmp_2_3 + tmp_4_3
        tmp_out_2_radix3 = tmp_2_3
        if (s1 == 0):
            self.output_2 = tmp_out_2_radix3
        elif (s1 == 1):
            self.output_2 = tmp_out_2_radix5
        
        self.output_4 = tmp_1_3 - tmp_3_3
        self.output_3 = tmp_2_3 - tmp_4_3

# radix 3 rotator
class MixedRadix_Rotator:
    cnt = 0
    def __init__(self, stage_index, size):
        self.input = 0.0
        self.output = 0.0
        self.stage_index = stage_index
        # self.twiddleROM = (np.ones(3**(num_of_stages-stage_index))).astype(complex)
        self.twiddleROM = (np.ones(size)).astype(complex)
        # self.two_thirds_len = 3**(num_of_stages-stage_index) - (3**(num_of_stages-stage_index)//3)
        self.two_thirds_len = size - size//3
        # print(self.two_thirds_len)
        # N = 3**num_of_stages
        N = size
        if (self.two_thirds_len > 2):
            for i in range(self.two_thirds_len):
                if (i < self.two_thirds_len//2):
                    k = i * 3**(stage_index)
                else:
                    k = 2*(i-self.two_thirds_len//2) * 3**(stage_index)
                # print("k = ", k)
                self.twiddleROM[i+(len(self.twiddleROM) - self.two_thirds_len)] = np.exp(-1j*2*np.pi*k/N)
        
        print(f"STAGE RADIX3, twiddle = {self.twiddleROM}")

        # print(self.twiddleROM)

    def rotate(self, fifo_full_flag):
        if (fifo_full_flag): # dozvola da brojac vrti i cita redom twiddle faktore iz memorije
            # print(f"STAGE {self.stage_index}, FIFO FULL")
            self.output = self.input * self.twiddleROM[self.cnt]
            # print(f'STAGE {self.stage_index}, curr twiddle = {self.twiddleROM[self.cnt]}')
            if (self.cnt == len(self.twiddleROM)-1):
                self.cnt = 0
            else:
                self.cnt += 1


class Radix5PhaseController:
    """
    Hardware-style phase generator for radix-5 SDF control.
    No division is used at runtime.

    phase:     0..4
    delay_cnt: 0..(cfg_delay-1)

    Each phase is held for cfg_delay valid samples.
    """
    def __init__(self, cfg, cfg_delay):
        self.cfg = int(cfg)
        self.cfg_delay = int(cfg_delay)
        if self.cfg_delay < 1:
            raise ValueError("cfg_delay must be >= 1")
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


class MixedRadix_SDF_stage_counter_ctrl:
    """
    Radix-2-3-5 SDF stage with hardware-like control:
    - input/output muxing based on phase counter
    - phase generated by (phase, delay_cnt) counters
    - cfg register supplies delay (no divide in datapath)

    Expected external classes already defined in this notebook:
      Fifo, MixedRadix_PreAdder, MixedRadix_Rotator
    """
    def __init__(self, config, stage_index, size, cfg_delay=None, s0=2, s1=1):
        if int(config) not in [2, 3, 5]:
            raise ValueError("This implementation is for radix-2-3-5 stages (config=5 or 3 or 2)")

        self.config = int(config)
        self.stage_index = int(stage_index)
        self.num_of_samples = int(size)

        # Config-register value for delay (programmed by software/control plane)
        # If not provided, default to size//5 for convenience.
        self.cfg_delay = int(cfg_delay) if cfg_delay is not None else (self.num_of_samples // 5)
        if self.cfg_delay < 1:
            raise ValueError("cfg_delay must be >= 1")

        # Optional preadder mode selects, kept configurable
        self.s0 = int(s0)
        self.s1 = int(s1)

        self.fifo_0 = Fifo(self.cfg_delay)
        self.fifo_1 = Fifo(self.cfg_delay)
        self.fifo_2 = Fifo(self.cfg_delay)
        self.fifo_3 = Fifo(self.cfg_delay)

        self.pre_adder = MixedRadix_PreAdder()
        self.rotator = MixedRadix_Rotator(stage_index=self.stage_index, size=self.num_of_samples)
        self.rot_en = 0

        self.ctrl = Radix5PhaseController(cfg=self.config, cfg_delay=self.cfg_delay)

        self.input_sample = 0.0
        self.output_sample = 0.0
        self.op_cnt = 0

    def set_muxes(self, s0, s1):
        self.s0 = int(s0)
        self.s1 = int(s1)

    def reset(self):
        self.ctrl.reset()
        self.input_sample = 0.0
        self.output_sample = 0.0
        self.op_cnt = 0

    def _load_preadder_inputs(self):
        self.pre_adder.input_0 = self.fifo_0.get_output()
        self.pre_adder.input_1 = self.fifo_1.get_output() if self.config != 2 else self.input_sample

        if self.config == 5:
            self.pre_adder.input_2 = self.fifo_2.get_output()
        elif self.config == 3:
            self.pre_adder.input_2 = self.input_sample
        else:  # self.config == 2
            self.pre_adder.input_2 = complex(0.0, 0.0)  # Not used in radix-2 mode

        self.pre_adder.input_3 = self.fifo_3.get_output() if self.config == 5 else complex(0.0, 0.0)
        self.pre_adder.input_4 = self.input_sample if self.config == 5 else complex(0.0, 0.0)

    def calculate(self, input_sample, valid=True):
        if not valid:
            return self.output_sample

        self.input_sample = input_sample

        # Combinational preadder path (always evaluated from current fifo outputs + new input)
        self._load_preadder_inputs()
        match self.config:
            case 2:
                self.s0 = 0
            case 3:
                self.s0 = 1
                self.s1 = 0
            case 5:
                self.s0 = 2
                self.s1 = 1
        self.pre_adder.calculate(s0=self.s0, s1=self.s1)

        # Muxing based only on current phase (no division)
        phase = self.ctrl.phase
        # print(f"Phase: {phase}, op_cnt: {self.op_cnt}")

        if phase == 0:
            pre_rot = self.fifo_0.get_output()
            self.fifo_0.shift(self.input_sample)
        elif phase == 1:
            pre_rot = self.fifo_1.get_output() if (self.config != 2) else self.pre_adder.output_0
            if (self.config != 2):
                self.fifo_1.shift(self.input_sample)
            else:
                self.fifo_0.shift(self.pre_adder.output_1)
        elif phase == 2:
            pre_rot = self.fifo_2.get_output() if (self.config != 3) else self.pre_adder.output_0
            if (self.config != 3):
                self.fifo_2.shift(self.input_sample)
            else:
                self.fifo_0.shift(self.pre_adder.output_1)
                self.fifo_1.shift(self.pre_adder.output_2)
        elif phase == 3:
            pre_rot = self.fifo_3.get_output()
            self.fifo_3.shift(self.input_sample)
        else:  # phase == 4
            pre_rot = self.pre_adder.output_0
            self.fifo_0.shift(self.pre_adder.output_1)
            self.fifo_1.shift(self.pre_adder.output_2)
            self.fifo_2.shift(self.pre_adder.output_3)
            self.fifo_3.shift(self.pre_adder.output_4)

        # Twiddle/rotator stage
        self.rotator.input = pre_rot
        match self.config:
            case 2:
                rot_en = self.fifo_0.is_full()  # Rotate when fifo_0 is full (phase 1)
            case 3:
                rot_en = self.fifo_1.is_full()  # Rotate when fifo_1 is full (phase 2)
            case 5:
                rot_en = self.fifo_3.is_full()  # Rotate when fifo_3 is full (phase 4)
        self.rotator.rotate(rot_en)
        self.output_sample = self.rotator.output

        # Advance counters (hardware tick)
        self.ctrl.tick(valid=True)

        if self.op_cnt == (self.num_of_samples - 1):
            self.op_cnt = 0
        else:
            self.op_cnt += 1

        return self.output_sample