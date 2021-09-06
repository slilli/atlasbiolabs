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
import hashing_ABO_permutations as h
import output as out
import filter_relevant_snps as snp
import matplotlib.pyplot as plt
# ----------------------------------------------------------

'''
@read_me:
    - hashing_ABO_permutations.py: dictionary for ABO permutations (A/B/AB/0) and genotypes (homozygous alternate/reference or heterozygous)
    - output.py: script to create .txt output, including samples, genotypes and bloodtypes
    - filter_relevant_snps.py: @here: 'rs8176719', 'rs8176746', 'rs8176747'
	- .vcf must be in the same directory as this python script
	- output file is created in the same directory as this python script
    - we assume,that snps in vcf files are on the forward strand, otherwise conversion is necessary
'''


def main(id_vcf, gt_vcf, samples_vcf, output_path):


    """
    This function filters out snps, computes each samples genotype associated bloodtype
    Finally, an output file is created containing interrelated sample_id, genotypes for all 3 relevant snps and the resulting blood type
    """


    # filter relevant snps
    x = snp.Snps(id_vcf, gt_vcf, samples_vcf, 'rs8176719', 'rs8176746', 'rs8176747')
    x.filter_snps()


    # compute genotye for all samples
    dict_gt = h.is_het_or_hom_dict()
    h.is_het_or_hom(dict_gt, x.gt_iterator, (np.arange(0, len(x.gt_subset[0]),1)), x.gt_hom_alt, x.gt_hom_ref, x.gt_het, x.gt_0, x.gt_1, x.gt_2)
    
    
    # computing blood types from genotype
    dict_obj = h.dictionary()
    h.hash_abo_permutations(dict_obj, x.blood_group, x.gt_0, x.gt_1, x.gt_2)


    # cretate txt output
    out.create_output(samples_vcf, x.gt_0, x.gt_1, x.gt_2, x.blood_group, output_path)
   

if __name__ == "__main__":


	# check if input criteria are met--------------
	try:
		if not len(sys.argv) == 3:
			raise ValueError("error")
	except(ValueError, IndexError):
		exit("\033[1m" + "2 Arguments needed: Your input should look like: [(...).py] [input_vcf] [output.file]")


	try:
		with open(sys.argv[1], mode='r') as vcf:
			print("opened vcf")
	except(IndexError):
		exit("Couldn't open vcf")


	# @param from .vcf----------------------
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']
	gt_from_vcf = callset['calldata/GT']
	samples_from_vcf = callset['samples']
	output_path = sys.argv[2]

	# calling functions ----------------------------
	main(id_from_vcf, gt_from_vcf, samples_from_vcf, output_path)
    

