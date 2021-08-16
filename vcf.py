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


# function to filter out relevant snps --------------------
def filtered_id(id_vcf, gt_vcf, samples_vcf, ref_vcf, alt_vcf):


	# check for relevant snps -------------------------
	values = ['rs8176719', 'rs8176746', 'rs8176747']
	boolean_array = np.in1d(values, id_vcf)
	print("Check for snps 'rs8176719', 'rs8176746', 'rs8176747:", boolean_array)

	
	try:
		for k in [0, len(boolean_array)-1]:
			if boolean_array[k]  == False :
				raise ValueError ("one(or more) required snps is not in variants/ID")
	except(ValueError, IndexError):
		exit("\033[1m" + "Could not determine ABO blood type. Not all required snps are present" + "\033[0m")
	

	# find indices of relevant snps in numpy array -------	
	value_index_19 = list(id_vcf).index('rs8176719')
	value_index_46 = list(id_vcf).index('rs8176746')
	value_index_47 = list(id_vcf).index('rs8176747')


	value_indices = [value_index_19, value_index_46, value_index_47]
	print("Position of relevant snps: ", value_indices)

	
	# create subset of genotype array (only including the three relevant snps)----
	gt_array = allel.GenotypeArray(gt_vcf)


	print("------genotype array subset------------------------------------------------")
	k = np.arange(0, len(samples_vcf), 1)
	gt_subset = gt_array.subset(value_indices, k)
	print(gt_subset)
	print("---------------------------------------------------------------------------")
	
	'''
	print("REF rs8176719: ", ref_vcf[value_index_19], "ALT rs8176719: ", alt_vcf[value_index_19])
	print("REF rs8176746: ", ref_vcf[value_index_46], "ALT rs8176746: ", alt_vcf[value_index_46])
	print("REF rs8176747: ", ref_vcf[value_index_47], "ALT rs8176747: ", alt_vcf[value_index_47])
	'''

		
	all_gt = []
	
	for i in [0, len(gt_subset)-1]:
		gt_hom_alt = gt_subset[i].is_hom_alt()
		gt_hom_ref = gt_subset[i].is_hom_ref()
		gt_het = gt_subset[i].is_het()
		for k in [0, len(gt_subset[i])-1]:
			if gt_hom_alt[k] == True:
				all_gt = np.append(all_gt, 'a')
			elif gt_hom_ref[k] == True:
				all_gt = np.append(all_gt, 'b')
			else:
				all_gt = np.append(all_gt, 'c')	

	
	print(all_gt)
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
	
	
	# @param --------------------------------------
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']
	gt_from_vcf = callset['calldata/GT']
	ref_from_vcf = callset['variants/REF']
	alt_from_vcf = callset['variants/ALT']
	

	# calling functions ----------------------------
	samples_from_vcf = callset['samples']


	filtered_id(id_from_vcf, gt_from_vcf, samples_from_vcf, ref_from_vcf, alt_from_vcf)


