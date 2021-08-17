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
import csv
# ----------------------------------------------------------

'''
@read_me:
	.vcf must be in the same directory as this python script
	.txt is created in the same directory as this python script
'''

# function to filter out relevant snps --------------------
def filtered_id(id_vcf, gt_vcf, samples_vcf, output_path):


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


	# find indices of relevant snps in numpy array and store them in value_indices-------
	value_index_19 = list(id_vcf).index('rs8176719')
	value_index_46 = list(id_vcf).index('rs8176746')
	value_index_47 = list(id_vcf).index('rs8176747')
	value_indices = [value_index_19, value_index_46, value_index_47]


	# create subset of genotype array (only including the three relevant snps)----
	gt_array = allel.GenotypeArray(gt_vcf)
	k = np.arange(0, len(samples_vcf), 1)
	gt_subset = gt_array.subset(value_indices, k)


	gt_hom_alt = gt_subset.is_hom_alt()
	gt_hom_ref = gt_subset.is_hom_ref()
	gt_het = gt_subset.is_het()
	gt_iterator = np.arange(0, len(gt_subset), 1)


	gt_19 = []
	gt_46 = []
	gt_47 = []


	for i in gt_iterator:
		for k in (np.arange(0, len(gt_subset[0]),1)):
			if gt_hom_alt[i][k] == True:
				if (i == 0):
					gt_19.append('rs8176719(I;I)')
				elif (i == 1):
					gt_46.append('rs8176746(C;C)')
				else:
					gt_47.append('rs8176747(G;G)')
			elif gt_hom_ref[i][k] == True:
				if (i == 0):
					gt_19.append('rs8176719(D;D)')
				elif (i == 1):
					gt_46.append('rs8176746(A;A)')
				else:
					gt_47.append('rs8176747(C;C)')
			else:
				if (i == 0):
					gt_19.append('rs8176719(D;I)')
				elif (i == 1):
					gt_46.append('rs8176746(A;C)')
				else:
					gt_47.append('rs8176747(C;G)')


	blood_group = []

	for it in (np.arange(0, len(gt_19), 1)):
		if ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(A;A)') and (gt_47[it] == 'rs8176747(C;C)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(A;A)') and (gt_47[it] == 'rs8176747(C;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] ==	'rs8176746(A;A)') and (gt_47[it] == 'rs8176747(G;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(A;C)') and (gt_47[it] == 'rs8176747(C;C)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(A;C)') and (gt_47[it] == 'rs8176747(C;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(A;C)') and (gt_47[it] == 'rs8176747(G;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(C;C)') and (gt_47[it] == 'rs8176747(C;C)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(C;C)') and (gt_47[it] == 'rs8176747(C;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;D)') and (gt_46[it] == 'rs8176746(C;C)') and (gt_47[it] == 'rs8176747(G;G)')):
			blood_group.append('Type 0')
		elif ((gt_19[it] == 'rs8176719(D;I)') and (gt_46[it] == 'rs8176746(C;C)') and (gt_47[it] == 'rs8176747(G;G)')):
			blood_group.append('Type A')
		elif ((gt_19[it] == 'rs8176719(I;I)') and (gt_46[it] == 'rs8176746(C;C)') and (gt_47[it] == 'rs8176747(G;G)')):
			blood_group.append('Type A')
		elif ((gt_19[it] == 'rs8176719(D;I)') and (gt_46[it] == 'rs8176746(A;A)') and (gt_47[it] == 'rs8176747(C;C)')):
			blood_group.append('Type B')
		elif ((gt_19[it] == 'rs8176719(I;I)') and (gt_46[it] == 'rs8176746(A;A)') and (gt_47[it] == 'rs8176747(C;C)')):
			blood_group.append('Type B')
		elif ((gt_19[it] == 'rs8176719(I;I)') and (gt_46[it] == 'rs8176746(A;C)') and (gt_47[it] == 'rs8176747(C;G)')):
			blood_group.append('Type AB')
		elif ((gt_19[it] == 'rs8176719(D;I)') and (gt_46[it] == 'rs8176746(A;C)') and (gt_47[it] == 'rs8176747(C;G)')):
			blood_group.append('Type B likely, Type A possible')
		else:
			blood_group.append('not clear')


# cretate txt output -----------------------------------------------------------


	with open(output_path,'w') as table:
		for row in zip(samples_vcf, gt_19, gt_46, gt_47, blood_group):
			for cell in row:
				table.write(str(cell) + '\t')
			table.write('\n')

	print(".txt containing blood types of all samples is in the same directory as your vcf")



if __name__ == "__main__":


	# check if input criteria are met--------------
	try:
		if not len(sys.argv) == 3:
			raise ValueError("error")
	except(ValueError, IndexError):
		exit("\033[1m" + "2 Arguments needed: Your input should look like: [(...).py] [input_vcf] [output.txt]")


	try:
		with open(sys.argv[1], mode='r') as vcf:
			print("opened vcf")
	except(IndexError):
		exit("Couldn't open vcf")


	# @param --------------------------------------
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']
	gt_from_vcf = callset['calldata/GT']
	samples_from_vcf = callset['samples']
	output_path = sys.argv[2]

	# calling functions ----------------------------
	filtered_id(id_from_vcf, gt_from_vcf, samples_from_vcf, output_path)
