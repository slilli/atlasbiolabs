import allel
import numpy as np

'''
@read_me:
    Class for storing attributes
'''


class Snps:

    """
    Initialise attributes that are used throughout the main.py script
    """

    def __init__(self, id_vcf, gt_vcf, samples_vcf, value0, value1, value2):
        
        #relevant snps and boolean_array to check if they are present
        self.values = [value0, value1, value2]
        self.boolean_array = np.in1d(self.values, id_vcf)
        
        #position(index) of the relevant snps in variants/ID
        self.value_index_0 = list(id_vcf).index(value0)
        self.value_index_1 = list(id_vcf).index(value1)
        self.value_index_2 = list(id_vcf).index(value2)
        self.value_indices = [self.value_index_0, self.value_index_1, self.value_index_2]
        
        #GenotypeArray with genotypes from all samples
        self.gt_array = allel.GenotypeArray(gt_vcf)
        
        #genotye array subset with genoptypes for relevant snps
        self.gt_subset = self.gt_array.subset(self.value_indices, np.arange(0, len(samples_vcf), 1))

        #heterozygous/homozygous        
        self.gt_hom_alt = self.gt_subset.is_hom_alt()
        self.gt_hom_ref = self.gt_subset.is_hom_ref()
        self.gt_het = self.gt_subset.is_het()
        
        #iterator 
        self.gt_iterator = np.arange(0, len(self.gt_subset), 1)
        
        #empty arrays to fill with genotype and blood group
        self.gt_0 = []
        self.gt_1 = []
        self.gt_2 = []
        self.blood_group = []


    """
    Check if all required snps are present in the .vcf
    """
    
    def filter_snps(self):
    
        print("Check for snps 'rs8176719', 'rs8176746', 'rs8176747':", self.boolean_array)
    
    
        try:
            for k in [0, len(self.boolean_array)-1]:
                if self.boolean_array[k]  == False :
                    raise ValueError ("one(or more) required snps is not in variants/ID")
        except(ValueError, IndexError):
            exit("\033[1m" + "Could not determine ABO blood type. Not all required snps are present" + "\033[0m")
