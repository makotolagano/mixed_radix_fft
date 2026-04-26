import numpy as np
from scipy.fft import fft
from mixed_radix_fft import *

# # Example: stage0 of radix-5 chain (N=25 -> delay=5)
# stage0_hw = MixedRadix_SDF_stage_counter_ctrl(config=5, stage_index=0, size=25, cfg_delay=5)

# # Example: stage1 of radix-5 chain (N=5 -> delay=1), works for stage_index != 0
# stage1_hw = MixedRadix_SDF_stage_counter_ctrl(config=5, stage_index=1, size=5, cfg_delay=1)

# # Quick sanity run
# x = np.append(np.arange(25), np.zeros(25)) 
# y = []
# for sample in x:
#     y0 = stage0_hw.calculate(sample, valid=True)
#     y1 = stage1_hw.calculate(y0, valid=True)
#     y.append(y1)

# print("Generated samples:", len(y))
# print("First 10 outputs:", np.array(y[:50]))

# # numpy fft
# fft_numpy = fft(x[:25])  # Only the first 25 samples are the actual input signal
# print ("Numpy FFT:", fft_numpy)

# Example: stage0 of radix-3 chain (N=9 -> delay=5)
stage0_hw = MixedRadix_SDF_stage_counter_ctrl(config=3, stage_index=0, size=9, cfg_delay=3)

# Example: stage1 of radix-3 chain (N=5 -> delay=1), works for stage_index != 0
stage1_hw = MixedRadix_SDF_stage_counter_ctrl(config=3, stage_index=0, size=3, cfg_delay=1)

# Quick sanity run
x = np.append(np.arange(9), np.zeros(10)) 
y = []
for sample in x:
    y0 = stage0_hw.calculate(sample, valid=True)
    y1 = stage1_hw.calculate(y0, valid=True)
    y.append(y1)

print("Generated samples:", len(y))
print("First 10 outputs:\n", np.array(y))

# numpy fft
fft_numpy = fft(x[:9])  # Only the first 9 samples are the actual input signal
print ("Numpy FFT:\n", fft_numpy)