#vcf file handling

#import modules ------------------

import allel
import sys
import numpy as np

#---------------------------------

def setup():


	if not len(sys.argv) == 2:
		print("Error: 2 arguments needed")
		print("[(...).py][input_vcf]")
		sys.exit()
	

	with open(sys.argv[1], mode='r') as vcf:
		print("opened vcf")


def extract_snps(s):

	
	print("Loading all snps into an array...")
	print(s)
	print("vcf_ID_length:", len(s))
	

def select_id():

	
	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = np.array(callset['variants/ID'])


	vcf_ID_new = []


	for k in id_from_vcf:
		#if (k[0] == 'r') and (k[1] == 's'):
			vcf_ID_new = np.append(id_from_vcf[k])
			print("vcf_ID_rs:", len(vcf_ID_new)) 


#main-----------------------------

if __name__ == "__main__":

	callset = allel.read_vcf(sys.argv[1])
	id_from_vcf = callset['variants/ID']
	setup()
	extract_snps(id_from_vcf)
	select_id()
