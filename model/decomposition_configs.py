import numpy as np

i_arr = []
j_arr = []
k_arr = []
for i in range(10):
    for j in range(10):
        for k in range(10):
            if ((((2**i) * (3**j) * (5**k)) <= 275)):
                i_arr.append(i)
                j_arr.append(j)
                k_arr.append(k)

N_arr = []
i_arr = np.unique(i_arr)
j_arr = np.unique(j_arr)
k_arr = np.unique(k_arr)
stage_nums = []

for i in i_arr:
    for j in j_arr:
        for k in k_arr:
            if((12 * ((2**i) * (3**j) * (5**k))) <= 3300):
            # if (not ((12 * (2**i * 3**j * 5**k)) in N_arr)):
                # print(f"i = {i}, j = {j}, k = {k}")
                N_arr.append(12 * (2**i * 3**j * 5**k))
                stage_nums.append(i+2 + j+1 + k)
                print("2: ", i+2, " | 3: ", j+1, " | 5: ", k)
                print(12 * (2**i * 3**j * 5**k))

N_arr = np.array(N_arr)
N_arr.sort()

print("Largest mixed-radix FFT sequence length = ", max(N_arr))
print("Number of radix combinations = ", len(N_arr))

print("Sequence lengths = ", N_arr)
print("Number of radix-2 stages = ", i_arr+2)
print("Number of radix-3 stages = ", j_arr+1)
print("Number of radix-5 stages = ", k_arr)
print("stage nums = ", np.array(stage_nums))