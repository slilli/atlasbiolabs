# ----------------------------------------------------------
# Python script to determine ABO blood group from a vcf file
# ATLAS Biolabs
# EPI-Dx Project
# Author: Lilli Schuckert
# ----------------------------------------------------------


#import modules --------------------------------------------
import allel
import sys
import numpy as np
import pandas as pd
# ----------------------------------------------------------



# check if input criteria are met ----------------------
def setup():


	if not len(sys.argv) == 2:
		print("Error: 2 arguments needed")
		print("Your input should look like: [(...).py][input_vcf]")
		sys.exit()


	print("--------------------------------------------------------------------------")


	with open(sys.argv[1], mode='r') as vcf:
		print("opened vcf")


# loading all snps into a numpy array -------------------
def extract_snps():


	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']


	print("--------------------------------------------------------------------------")


	print("Loading all snps into an array...")
	print(id_from_vcf)
	print("vcf_ID_length:", len(id_from_vcf))


# function to filter out relevant snps --------------------
def filtered_id():


	# variants/ID -------------------------------------
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']


	print("--------------------------------------------------------------------------")
	

	# check for relevant snps -------------------------
	values = ['rs8176747' , 'rs8176746', 'rs8176719']
	boolean_array = np.in1d(values, id_from_vcf)
	print("Check for snps 'rs8176747', 'rs8176746', 'rs8176719':", boolean_array)


	try:
		for k in [0, len(boolean_array)]:
			if boolean_array[k]  == False :
				raise ValueError ("one(or more) required snps is not in variants/ID")
	except(ValueError, IndexError):
		exit("\033[1m" + "Could not determine ABO blood type. Not all required snps are present" + "\033[0m")
	
	'''	
	result_1 = np.where(id_from_vcf == 'rs8176747')
	print("snp rs8176747 is located at position" ,result_1, "in the numpy array")	

	result_2 = np.where(id_from_vcf == 'rs8176746')
	print("snp rs8176746 is located at position" ,result_2, "in the numpy array")
	
	result_3 = np.where(id_from_vcf == 'rs8176719')
	print("snp rs8176719 is located at position" ,result_3, "in the numpy array")
	'''

#main-----------------------------


if __name__ == "__main__":


	setup()
	extract_snps()
	filtered_id()


