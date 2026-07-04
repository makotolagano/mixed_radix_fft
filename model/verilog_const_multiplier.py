#!/usr/bin/python
# Author: Karthikeyan
# https://github.com/lxkarthi
# Description: Generates Verilog code for fractional constant multiplier
## Functions defined
##############################################################################################
def get_value_prompt(prompt, vartype, condition, condition_value, error_string, default_value):
	try:
		var = vartype(input(prompt))
		if not(condition(var,condition_value)):
			raise ValueError()
	except ValueError:
		print (error_string)
		print ("default=",default_value)
		var=default_value
	else:
		var=var
	return var
##############################################################################################

#FIXME may be not needed
# Total number of bits to represent the number
nobs_multiplier=32 
try:
	nob = int(input("Number of bits to represent:"))
	if (nob<1):
		raise ValueError()
except ValueError:
	print ("Invalid number of bits (should be >0)...")
	print ("default=",nobs_multiplier)
else:
	nobs_multiplier=nob

# Number of bits after decimal point
nobs_after_dp=nobs_multiplier 
try:
	nob = int(input("Number of bits after decimal point:"))
	if (nob>nobs_multiplier):
		raise ValueError()
except ValueError:
	print ("Invalid number of bits after decimal point(should be <",nobs_multiplier,")...")
	print ("default=",nobs_multiplier)
else:
	nobs_after_dp=nob

nobs_before_dp=nobs_multiplier-nobs_after_dp
integer_bits=nobs_multiplier-nobs_after_dp

# Multiplier
print (nobs_multiplier,"=",integer_bits,".",nobs_after_dp, ' '.join(map(str, ["Range(",2**integer_bits,"< value <",2**(-nobs_after_dp),")"])))
multiplier=get_value_prompt("Enter the number to multiply:", float, (lambda a, b: a < 2**integer_bits and a>=2**(-nobs_after_dp)), 0, ' '.join(map(str, ["Number out of range(",2**integer_bits,"< value <",2**(-nobs_after_dp),")"])), 0);
integer_part=int(multiplier)
fraction_part=multiplier-integer_part
actual_fraction=int(fraction_part*(2**nobs_after_dp))/float(2**nobs_after_dp)
print ("Actual multiplier=",integer_part+actual_fraction)
# Output multiplicand number of bits
nobs_multiplicand=get_value_prompt("Enter the number of bits in input multiplicand:", int, (lambda a, b: a > b), 0, "Invalid  number of bits(should be >0)...", nobs_multiplier);
print ("input[",nobs_multiplicand,"]*",integer_part+actual_fraction,"[",nobs_after_dp,".",nobs_after_dp,"]=output[",nobs_before_dp+nobs_multiplicand,".",nobs_after_dp,"]\n")

#Choice of multiplier
choice=get_value_prompt("Multiplier type\n1.*\n2.>>>\nChoose:", int, (lambda a, b:1 <= a <= 2), 0, "Invalid  choice...", 1);
print ("/////// VERILOG CODE FOR MULTIPLIER input[",nobs_multiplicand,"]*const[",nobs_after_dp,".",nobs_after_dp,"]=out[",nobs_before_dp+nobs_multiplicand,".",nobs_after_dp,"] /////////")
# Generate Verilog code

print ("module unsigned_mult (out, a, b);")
print ("	parameter outwidth = ",nobs_multiplier+nobs_multiplicand,", inwidth = ",nobs_multiplicand,";")
print ("	output [outwidth-1:0] out;")
print ("	input  [inwidth-1:0] a;")
# version 1
if(choice==1):
	print ("	assign out = a * ",int(multiplier*(2**nobs_after_dp)),";")
# version 2
elif(choice==2):
	bin_len = len(bin(int(multiplier*(2**nobs_after_dp))).replace("0b","")[::-1])
	print ("	assign out = ",end="")
	#print bin(int(multiplier*(2**nobs_after_dp))).replace("0b","")
	for (i,bit) in enumerate(bin(int(multiplier*(2**nobs_after_dp))).replace("0b","")[::-1]):
		if(bit=='1'):
			print ("a<<<",i,end=" ")
			if(i!= bin_len-1):
				print (" + ",end=" ")
	print (";")
print ("endmodule")
