import numpy as np

'''
@read_me:
   - two dictionaries: ABO blood group permutations and a dictionary for genotypes as well as two functions to access
     those dictionaries
'''


def dictionary():

    """
    Nested dictionary with key-{key1:value1,....,key4-value4} pairs, where if the first 3 values match the 3 snps from our vcf file, the 4th value
    corresponds to the resulting bloodtype
    """
    
    f = open('ABO_bloodgroup_permutations.txt')

    # Create empty dictionary 
    abo_bloodtype_dict = {}

    for line in f:
        x = line.rstrip().split("  ")
        if not x[0] in abo_bloodtype_dict:
            abo_bloodtype_dict[x[0]] = {x[1] : {x[2] : x[3]}}
        else:
            if not x[1] in abo_bloodtype_dict[x[0]]:
                abo_bloodtype_dict[x[0]][x[1]] = {}
                abo_bloodtype_dict[x[0]][x[1]][x[2]] = x[3]
            else:
                abo_bloodtype_dict[x[0]][x[1]][x[2]] = x[3]
    
    return abo_bloodtype_dict

'''
def dictionary():

    """
    Nested dictionary with key-{key1:value1,....,key4-value4} pairs, where if the first 3 values match the 3 snps from our vcf file, the 4th value
    corresponds to the resulting bloodtype
    """
    
    abo_bloodtype_dict = {
    'rs8176719(D;D)':{
        'rs8176746(A;A)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'Type O'},
            'rs8176747(C;G)':{
                'Blood Type' : 'Type O'},
            'rs8176747(G;G)':{
                'Blood Type' : 'Type O'}},
        'rs8176746(A;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'Type O'},
            'rs8176747(C;G)':{
                'Blood Type' : 'Type O'},
            'rs8176747(G;G)':{
                'Blood Type' : 'Type O'}},
        'rs8176746(C;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'Type O'},
            'rs8176747(C;G)':{
                'Blood Type' : 'Type O'},
            'rs8176747(G;G)':{
                'Blood Type' : 'Type O'}}},
    'rs8176719(D;I)':{
        'rs8176746(A;A)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'Type B'},
            'rs8176747(C;G)':{
                'Blood Type' : 'not clear'},
            'rs8176747(G;G)':{
                'Blood Type' : 'not clear'}},
        'rs8176746(A;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'not clear'},
            'rs8176747(C;G)':{
                'Blood Type' : 'Type B likely, Type A possible'},
            'rs8176747(G;G)':{
                'Blood Type' : 'not clear'}},
        'rs8176746(C;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'not clear'},
            'rs8176747(C;G)':{
                'Blood Type' : 'not clear'},
            'rs8176747(G;G)':{
                'Blood Type' : 'Type A'}}},
    'rs8176719(I;I)':{
        'rs8176746(A;A)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'Type B'},
            'rs8176747(C;G)':{
                'Blood Type' : 'not clear'},
            'rs8176747(G;G)':{
                'Blood Type' : 'not clear'}},
        'rs8176746(A;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'not clear'},
            'rs8176747(C;G)':{
                'Blood Type' : 'Type AB'},
            'rs8176747(G;G)':{
                'Blood Type' : 'not clear'}},
        'rs8176746(C;C)':{
            'rs8176747(C;C)':{
                'Blood Type' : 'not clear'},
            'rs8176747(C;G)':{
                'Blood Type' : 'not clear'},
            'rs8176747(G;G)':{
                'Blood Type' : 'Type A'}}}}
    
    
    return abo_bloodtype_dict
'''

def hash_abo_permutations(bloodtype_dict, blood_gr, rs_array0, rs_array1, rs_array2):

    """
    Check wether the genotype from our sampls are matching one of the key-value pairs from the dictionary
    """
    
    for k in (np.arange(0, len(rs_array0),1)):
        blood_gr.append(bloodtype_dict[rs_array0[k]][rs_array1[k]][rs_array2[k]])
            
    return blood_gr, rs_array0, rs_array1, rs_array2


def is_het_or_hom_dict():

    """
    Dictionary for heterozygus/homozygous genotype
    """

    het_hom_dict = {
    0:{
        'hom_alt':{'snp' : 'rs8176719(I;I)'},
        'hom_ref': {'snp' : 'rs8176719(D;D)'},
        'het' : {'snp' : 'rs8176719(D;I)'}},
    1:{
        'hom_alt':{'snp' : 'rs8176746(C;C)'},
        'hom_ref':{'snp' : 'rs8176746(A;A)'},
        'het':{'snp' : 'rs8176746(A;C)'}},
    2:{
        'hom_alt':{'snp' : 'rs8176747(D;D)'},
        'hom_ref':{'snp' : 'rs8176747(B;B)'},
        'het':{'snp' : 'rs8176747(B;D)'}}}

    
    return het_hom_dict


def is_het_or_hom(het_hom_dict, num_gt_subset, len_gt_subset, gt_hom_alt, gt_hom_ref, gt_het, rs_array0, rs_array1, rs_array2):
    
    """
    Check if genotype is heterozygous/homozygous alternate/homozygous reference. Access dictionary at corresponding position
    """

    arrays = [rs_array0, rs_array1, rs_array2]
    
    
    for i in num_gt_subset:
        for k in len_gt_subset:
            if (gt_hom_alt[i][k] == True):           
                arrays[i].append(het_hom_dict[i]['hom_alt']['snp'])
            elif (gt_hom_ref[i][k] == True):
                arrays[i].append(het_hom_dict[i]['hom_ref']['snp'])
            else:
                arrays[i].append(het_hom_dict[i]['het']['snp'])
                    

    return arrays[0], arrays[1], arrays[2]
