'''
@read_me:
    creates an output file by zipping all desired arrays together
    specify your file type when running main.py
'''

def create_output(samples, rs_array0, rs_array1, rs_array2, blood_group_array, output_path):

    with open(output_path,'w') as table:
        for row in zip(samples, rs_array0, rs_array1, rs_array2, blood_group_array):
            for cell in row:
                table.write(str(cell) + '\t')
            table.write('\n')

    print("The output file containing blood types of all samples is in the same directory as your .vcf")