import numpy as np
import fxpmath as fxp

class Fifo:
    full = 0
    cnt = 0

    def __init__(self, depth):
        self.depth = depth
        self.buffer = (np.zeros(depth)).astype(complex)
        # print(self.depth)

    def is_full(self):
        return self.full
    
    def get_output(self):
        return self.buffer[-1]
    
    def shift(self, input_sample):
        self.cnt += 1
        self.buffer = np.roll(self.buffer, 1)
        # print("FIFO input samples = ", input_sample)
        self.buffer[0] = input_sample
        if (self.cnt > self.depth):
            self.full = 1
        else:
            self.full = 0

class Fifo_fxp:
    full = 0
    cnt = 0

    def __init__(self, depth, DATA=fxp.Fxp(None, 16, 15)):
        self.depth = depth
        self.DATA  = DATA
        self.buffer = Fxp(np.zeros((depth,2))).like(DATA)
        # print(self.depth)

    def is_full(self):
        return self.full
    
    def get_output(self):
        return self.buffer[-1]
    
    def shift(self, input_sample_r, input_sample_i):
        self.cnt += 1
        # print("FIFO buffer before roll = ", self.buffer)
        self.buffer = np.roll(self.buffer, 1, axis=0).like(self.DATA)
        # print("FIFO input samples = ", input_sample_r, input_sample_i)
        self.buffer[0,0] = input_sample_r
        self.buffer[0,1] = input_sample_i
        # print("FIFO buffer after roll = ", self.buffer)
        if (self.cnt > self.depth):
            self.full = 1
        else:
            self.full = 0

def digit_reverse(radices):
    radices_rev = np.copy(radices)
    radices_rev = radices_rev[::-1]

    mr_fft_len = np.prod(radices)
    indices = np.arange(mr_fft_len)

    mult_factors = []
    tmp_len = mr_fft_len
    for radix in radices:
        tmp_len = tmp_len//radix
        mult_factors.append(tmp_len)
    
    mult_factors_rev = []
    tmp_len = mr_fft_len
    for radix in radices_rev:
        tmp_len = tmp_len//radix
        mult_factors_rev.append(tmp_len)
    

    decomp = []
    for index in indices:
        tmp_decomposition = []
        tmp_index = index
        for mult_factor in mult_factors:
            tmp_decomposition.append(tmp_index // mult_factor)
            tmp_index = tmp_index % mult_factor
        decomp.append(tmp_decomposition)
    
    decomp = np.array(decomp)
    decomp = decomp.T[::-1]

    indices_rev = np.multiply(np.tile(mult_factors_rev, (mr_fft_len,1)), decomp.T).sum(axis=1)

    return indices_rev