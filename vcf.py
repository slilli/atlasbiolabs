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
	values = ['rs8176719', 'rs8176746', 'rs8176747']
	boolean_array = np.in1d(values, id_from_vcf)
	print("Check for snps 'rs8176719', 'rs8176746', 'rs8176747:", boolean_array)

	
	try:
		for k in [0, len(boolean_array)-1]:
			if boolean_array[k]  == False :
				raise ValueError ("one(or more) required snps is not in variants/ID")
	except(ValueError, IndexError):
		exit("\033[1m" + "Could not determine ABO blood type. Not all required snps are present" + "\033[0m")
	

	# find indices of relevant snps in numpy array -------	
	value_index_19 = list(id_from_vcf).index('rs8176719')
	print("snp rs8176719 is located at position: ", value_index_19)

	
	value_index_46 = list(id_from_vcf).index('rs8176746')
	print("snp rs8176746 is located at position: ", value_index_46)


	value_index_47 = list(id_from_vcf).index('rs8176747')
	print("snp rs8176747 is located at position: ", value_index_47)

	value_indices = [value_index_19, value_index_46, value_index_47]
	print(value_indices)

	
	# create genotype array from  all samples (including only the genotypes respectively to our snps)
	gt_array = allel.GenotypeArray(callset['calldata/GT'])

	print("length gt_array: ", len(gt_array))


	new = [gt_array[value_index_19], gt_array[value_index_46], gt_array[value_index_47]]
	print(new)


	ref_vcf = callset['variants/REF']
	alt_vcf = callset['variants/ALT']
	print("REF rs8176719: ", ref_vcf[value_index_19], "ALT rs8176719: ", alt_vcf[value_index_19])
	print("REF rs8176746: ", ref_vcf[value_index_46], "ALT rs8176746: ", alt_vcf[value_index_46])
	print("REF rs8176747: ", ref_vcf[value_index_47], "ALT rs8176747: ", alt_vcf[value_index_47])

	
if __name__ == "__main__":


	# check if input criteria are met--------------
	try:
		if not len(sys.argv) == 2:
			raise ValueError("error")
	except(ValueError, IndexError):
		exit("\033[1m" + "2 Arguments needed: Your input should look like: [(...).py] [input_vcf]")

	
	try:
		with open(sys.argv[1], mode='r') as vcf:
			print("opened vcf")
	except(IndexError):
		exit("Couldn't open vcf")
	
	
	# calling functions ---------------------------
	extract_snps()
	filtered_id()


