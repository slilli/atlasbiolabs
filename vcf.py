#vcf file handling
#import modules ------------------


import allel
import sys
import numpy as np
import pandas as pd


#---------------------------------


def setup():


	if not len(sys.argv) == 2:
		print("Error: 2 arguments needed")
		print("Your input should look like: [(...).py][input_vcf]")
		sys.exit()
	

	with open(sys.argv[1], mode='r') as vcf:
		print("opened vcf")


def extract_snps():


	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']


	print("Loading all snps into an array...")
	print(id_from_vcf)
	print("vcf_ID_length:", len(id_from_vcf))


def filtered_id():

	
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']
	ref_from_vcf = callset['variants/REF']
	alt_from_vcf = callset['variants/ALT']
	

	relevant_ID_for_ABO = []


	counter = 0
	counter2 = 0


	for k in [0,len(id_from_vcf)-1]:
		if(id_from_vcf[k] == "rs8176747") or (id_from_vcf[k] == "rs8176747") or (id_from_vcf[k] == "rs8176719"):
			counter2 += 1
			relevant_ID_for_ABO = np.append(relevant_ID_for_ABO, id_from_vcf[k])		


	print("relevant_ID_length:", len(relevant_ID_for_ABO)) 
	print(relevant_ID_for_ABO)
	print("counter:", counter)
	print("counter2:", counter2)


#main-----------------------------


if __name__ == "__main__":


	setup()
	extract_snps()
	filtered_id()


