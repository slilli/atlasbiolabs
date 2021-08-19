import numpy as np


def dictionary():

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
    '25': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;C)', 'snp3' : 'rs8176747(C;G)', 'Blood Type' : 'AB'},
    '26': {'snp1' : 'rs8176719(I;I)', 'snp2' : 'rs8176746(A;A)', 'snp3' : 'rs8176747(C;C)', 'Blood Type' : 'Type B'}}
    
    
    return abo_bloodtype_dict


def hash_abo_permutations(dict_obj, blood_gr, *rs_array):

    for k in (np.arange(0, len(rs_array[0]),1)):
        for i in dict_obj:
            if (rs_array[0][k] == dict_obj[i].get('snp1')) and (rs_array[1][k] == dict_obj[i].get('snp2')) and (rs_array[2][k] == dict_obj[i].get('snp3')) :
                blood_gr.append(dict_obj[i].get('Blood Type'))
    
    return blood_gr

                





























'''
datalists = dict()
for item in ['size', 'type']:
    datalists[item] = []


def nested_dict_values_iterator(dict_obj):

    # Iterate over all values of given dictionary
    for value in dict_obj.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over all its values
            for v in  nested_dict_values_iterator(value):
                #yield (value.get('Blood Type'))
                print(1)
                yield v
        else:
            # If value is not dict type then yield the value
            print(2)
            yield value

   
            
            
#Loop through all nested dictionary values
for value in nested_dict_values_iterator(abo_bloodtype_dict):
    print(value)
   


for d in abo_bloodtype_dict.values():
    if (d['snp1'] == 'rs8176719(I;I)') and (d['snp2'] == 'rs8176746(A;A)') and (d['snp3'] =='rs8176747(C;C)'):
        print("Yes")
        break
else:
    print("No")
'''