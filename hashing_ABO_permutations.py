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
    
    abo_bloodtype_dict = {
    '0': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type O'},
    '1': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'Type O'},
    '2': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'Type O'},
    '3': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type O'},
    '4': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'Type O'},
    '5': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'Type O'},
    '6': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type O'},
    '7': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'Type O'},
    '8': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'Type O'},
    '9': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'not clear'},
    '10': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'not clear'},
    '11': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'not clear'},
    '12': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'not clear'},
    '13': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'not clear'},
    '14': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'not clear'},
    '15': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'Type B likely, Type A possible'},
    '16': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'Type A'},
    '17': {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type B'},
    '18': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'not clear'},
    '19': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'not clear'},
    '20': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'not clear'},
    '21': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'not clear'},
    '22': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'not clear'},
    '23': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'not clear'},
    '24': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(G;G)', 'Blood Type' : 'Type A'},
    '25': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'Type AB'},
    '26': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type B'}}
    
    
    return abo_bloodtype_dict


def hash_abo_permutations(dict_obj, blood_gr, rs_array0, rs_array1, rs_array2):

    """
    Check wether the genotype from our sampls are matching one of the key-value pairs from the dictionary
    """
    for k in (np.arange(0, len(rs_array0),1)):
        for i in dict_obj:
            if (rs_array0[k] == dict_obj[i].get('snp1')) and (rs_array1[k] == dict_obj[i].get('snp2')) and (rs_array2[k] == dict_obj[i].get('snp3')) :
                blood_gr.append(dict_obj[i].get('Blood Type'))
            
    
    return blood_gr, rs_array0, rs_array1, rs_array2


def is_het_or_hom_dict():

    """
    Dictionary for heterozygus/homozygous genotype
    """

    het_hom_dict = {
    'hom_alt': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(C;C)', 'snp3' : 'rs8176747(G;G)'},
    'hom_ref': {'snp1' : 'rs8176719(D;D)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;C)'},
    'het' : {'snp1' : 'rs8176719(D;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;G)'}}
    
    return het_hom_dict


def is_het_or_hom(het_hom_dict, num_gt_subset, len_gt_subset, gt_hom_alt, gt_hom_ref, gt_het, rs_array0, rs_array1, rs_array2):
    
    """
    Check if genotype is heterozygous/homozygous alternate/homozygous reference. Access dictionary at corresponding position
    """
    
    for i in num_gt_subset:
        for k in len_gt_subset:
            if (gt_hom_alt[i][k] == True):
                if(i == 0):            
                    rs_array0.append(het_hom_dict['hom_alt'].get('snp1'))
                elif(i == 1):
                    rs_array1.append(het_hom_dict['hom_alt'].get('snp2'))
                else:
                    rs_array2.append(het_hom_dict['hom_alt'].get('snp3'))
            elif (gt_hom_ref[i][k] == True):
                if(i == 0):
                    rs_array0.append(het_hom_dict['hom_ref'].get('snp1'))
                elif(i == 1):
                    rs_array1.append(het_hom_dict['hom_ref'].get('snp2'))
                else:
                    rs_array2.append(het_hom_dict['hom_ref'].get('snp3'))
            else:
                if(i == 0):
                    rs_array0.append(het_hom_dict['het'].get('snp1'))
                elif(i == 1):
                    rs_array1.append(het_hom_dict['het'].get('snp2'))
                else:
                    rs_array2.append(het_hom_dict['het'].get('snp3'))
                    
    
    return rs_array0, rs_array1, rs_array2
