"""
Simple PreAdder Test - All 3 Radix Configurations
Tests with 3 simple input vectors and prints all intermediate stages
Shows both floating point and binary (fixed-point) representations
Uses same fixed-point format as the actual PreAdder model
"""

from mixed_radix_fft_fxp import MixedRadix_PreAdder_FXP
from fxpmath import Fxp

def get_binary_representation(val, dtype='fxp-s32/12'):
    """Get actual binary representation using Fxp library"""
    # Parse dtype to get bit width (e.g., 'fxp-s32/12' -> word_bits=32, frac_bits=12)
    import re
    match = re.match(r'fxp-s(\d+)/(\d+)', dtype)
    if match:
        word_bits = int(match.group(1))
        frac_bits = int(match.group(2))
    else:
        word_bits = 32
        frac_bits = 12
    
    if isinstance(val, complex):
        real_fxp = Fxp(val.real, True, dtype=dtype)
        imag_fxp = Fxp(val.imag, True, dtype=dtype)
        # Get the raw integer value in fixed-point format
        real_int = real_fxp.val
        imag_int = imag_fxp.val
        # Mask to word_bits and pad with zeros
        mask = (1 << word_bits) - 1
        real_bin = bin(real_int & mask)[2:].zfill(word_bits)
        imag_bin = bin(imag_int & mask)[2:].zfill(word_bits)
        return real_bin, imag_bin
    else:
        fxp_val = Fxp(val, True, dtype=dtype)
        fxp_int = fxp_val.val
        mask = (1 << word_bits) - 1
        binary = bin(fxp_int & mask)[2:].zfill(word_bits)
        return binary, None

def print_complex(val, dtype='fxp-s32/12'):
    """Format complex number with float and binary representation"""
    real_bin, imag_bin = get_binary_representation(val, dtype)
    
    float_str = f"{val.real:18.12f}+{val.imag:18.12f}j"
    binary_str = f"Real: {real_bin} | Imag: {imag_bin}"
    
    return f"{float_str:<22} [{binary_str}]"

def test_preadder():
    """Test PreAdder for all three radix configurations"""
    
    dtype = 'fxp-s17/12'
    
    # Three simple test input vectors
    test_inputs = [
        {
            'name': 'Zeros',
            'values': [0+0j, 0+0j, 0+0j, 0+0j, 0+0j],
        },
        {
            'name': 'Ones',
            'values': [1+0j, 1+0j, 1+0j, 1+0j, 1+0j],
        },
        {
            'name': 'Simple',
            'values': [1+1j, 0.5+0.5j, 0.3+0.3j, 0.2+0.2j, 0.1+0.1j],
        },
    ]
    
    configs = [
        {'radix': 2, 's0': 0, 's1': 0},
        {'radix': 3, 's0': 1, 's1': 0},
        {'radix': 5, 's0': 2, 's1': 1},
    ]
    
    for config in configs:
        print("\n" + "="*140)
        print(f"RADIX-{config['radix']} TEST (s0={config['s0']}, s1={config['s1']})")
        print("="*140)
        
        for test_case in test_inputs:
            print(f"\nInput: {test_case['name']}")
            print("-" * 140)
            
            preadder = MixedRadix_PreAdder_FXP(dtype=dtype, shift_bits=0)
            
            # Load inputs
            preadder.input_0 = test_case['values'][0]
            preadder.input_1 = test_case['values'][1]
            preadder.input_2 = test_case['values'][2]
            preadder.input_3 = test_case['values'][3]
            preadder.input_4 = test_case['values'][4]
            
            print(f"  input_0: {print_complex(preadder.input_0, dtype)}")
            print(f"  input_1: {print_complex(preadder.input_1, dtype)}")
            print(f"  input_2: {print_complex(preadder.input_2, dtype)}")
            print(f"  input_3: {print_complex(preadder.input_3, dtype)}")
            print(f"  input_4: {print_complex(preadder.input_4, dtype)}")
            
            # Calculate
            preadder.calculate(s0=config['s0'], s1=config['s1'])
            
            # Print outputs
            print(f"\n  OUTPUTS:")
            print(f"  output_0: {print_complex(preadder.output_0, dtype)}")
            print(f"  output_1: {print_complex(preadder.output_1, dtype)}")
            print(f"  output_2: {print_complex(preadder.output_2, dtype)}")
            print(f"  output_3: {print_complex(preadder.output_3, dtype)}")
            print(f"  output_4: {print_complex(preadder.output_4, dtype)}")

if __name__ == '__main__':
    test_preadder()
