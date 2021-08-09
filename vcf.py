## Lilli Schuckert
## Berufspraktikum Atlas Biolabs
## 03.08.2021
##------------------------------------------------------------------------------
import allel
import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
print(allel.__version__)
##------------------------------------------------------------------------------

##-------------- read vcf-------------------------------------------------------
with open('test.vcf',mode='r') as vcf:
    print(vcf.read())
##------------------------------------------------------------------------------

'''
The callset object returned by read_vcf() is a Python dictionary (dict).
It contains several NumPy arrays, each of which can be accessed via a key.
'''
callset = allel.read_vcf('test.vcf')
sorted(callset.keys())

callset['samples']
callset['variants/CHROM']
callset['variants/POS']
callset['variants/QUAL']

'''
All arrays with keys beginning ‘calldata/’ come from the sample fields in the
VCF file. For example, here are the actual genotype calls from the ‘GT’ field:
'''
callset['calldata/GT']


##-----------Aside: genotype arrays---------------------------------------------
gt = allel.GenotypeArray(callset['calldata/GT'])
print(gt)

'''
the is_het() method locates all heterozygous genotype calls:
'''
gt.is_het()

'''
the count_het() method will count heterozygous calls, summing over variants
(axis=0) or samples (axis=1) if requested
'''
gt.count_het(axis=1)

'''
here is how to perform an allele count, i.e., count the number times each allele
(0=reference, 1=first alternate, 2=second alternate, etc.) is observed for each variant:
'''
ac = gt.count_alleles()
print(ac)

##-------------Fields-----------------------------------------------------------
'''
VCF files can often contain many fields of data, and you may only need to extract
some of them to perform a particular analysis. You can select which fields to
extract by passing a list of strings as the fields parameter. For example,
let’s extract the ‘DP’ field from within the ‘INFO’ field, and let’s also extract
the ‘DP’ field from the genotype call data:
'''

callset = allel.read_vcf('test.vcf', fields=['variants/DP', 'calldata/DP'])
sorted(callset.keys())

##data we just extracted
callset['variants/DP']
callset['calldata/DP']

##extracting everything
callset = allel.read_vcf('test.vcf', fields='*')
sorted(callset.keys())

##----------Types---------------------------------------------------------------
##loaded into a 32-bit integer array
callset = allel.read_vcf('test.vcf', fields=['DP'])
callset['variants/DP']

##override to load data into 16-bit integer array
callset = allel.read_vcf('test.vcf', fields=['DP'], types={'DP': 'int16'})
callset['variants/DP']

##choose floating-point data type even for fields that are declares as type integer
callset = allel.read_vcf('test.vcf', fields=['DP'], types={'DP': 'float32'})
callset['variants/DP']

## usually scikit-allel will use an 'object' data type
'''
The advantage of using ‘object’ dtype is that strings can be of any length.
Alternatively, you can use a fixed-length string dtype, e.g.:
'''
callset = allel.read_vcf('test.vcf', types={'REF': 'S3'})
callset['variants/REF']
'''
Note that fixed-length string dtypes will cause any string values longer than
the requested number of characters to be truncated. I.e., there can be some data
loss. E.g., if using a single-character string for the ‘REF’ field, the correct
value of ‘GTC’ for the final variant will get truncated to ‘G’:
'''

##----------Numbers-------------------------------------------------------------
'''
for fields like ALT scikit uses default number of expected values which is set
at 3.
If you need to increase or decrease the expected number of values for any field,
you can do this via the numbers parameter. E.g., increase the number of ALT values to 5:
'''
#default
callset = allel.read_vcf('test.vcf')
callset['variants/ALT']
#increased
callset = allel.read_vcf('test.vcf', numbers={'ALT': 5})
callset['variants/ALT']

'''
Often there will be several fields within a VCF that all have a number of values
that depends on the number of alternate alleles (declared with number ‘A’ or ‘R’
in the VCF meta-information). You can set the expected number of values
simultaneously for all such fields via the alt_number parameter:
'''
callset = allel.read_vcf('test.vcf', fields=['ALT', 'AF'], alt_number=2)
callset['variants/ALT']
callset['variants/AF']

##---------Genotype ploidy------------------------------------------------------
'''
per default diploid organisms
'''
##-------Region-----------------------------------------------------------------
'''
You can extract data for only a specific chromosome or genome region via the
region parameter. The value of the parameter should be a region string of the
format ‘{chromosome}:{begin}-{end}’
'''
callset = allel.read_vcf('test.vcf', region='20:1000000-1231000')
callset['variants/POS']

##-------Samples----------------------------------------------------------------
'''
extract data for specific samples
'''
callset = allel.read_vcf('test.vcf', samples=['NA00001', 'NA00003'])
callset['samples']

allel.GenotypeArray(callset['calldata/GT']) # genotype array will only have 2 columns

## vcf_to_npz()-----------------------------------------------------------------
'''
NumPy arrays are stored in main memory (a.k.a., RAM), which means that as soon as
you end your Python session or restart your Jupyter notebook kernel, any data
stored in a NumPy array will be lost. If your VCF file is not too big, you can
extract data from the file into NumPy arrays then save those arrays to disk via
the vcf_to_npz() function. This function has most of the same parameters as the
read_vcf() function, except that you also specify an output path, which is the
name of the file you want to save the extracted data to.
'''

# only needs to be done once:
# allel.vcf_to_npz('example.vcf', 'example.npz', fields='*', overwrite=True)
# then:
# callset = np.load('example.npz')
# callset

## vcf_to_hdf5() for large datasets---------------------------------------------
'''
For large datasets, the vcf_to_hdf5() function is available. This function again
takes similar parameters to read_vcf(), but will store extracted data into an
HDF5 file stored on disk. The extraction process works through the VCF file in
chunks, and so the entire dataset is never loaded entirely into main memory
'''
# allel.vcf_to_hdf5('example.vcf', 'example.h5', fields='*', overwrite=True)
# import h5py
# callset = h5py.File('example.h5', mode='r')
# callset

'''
The one difference to be aware of here is that accessing data via a key like
‘variants/POS’ does not return a NumPy array, instead you get an HDF5
dataset object.
'''

## vcf_to_dataframe() ----------------------------------------------------------
df = allel.vcf_to_dataframe('test.vcf')
print(df)

df2 = allel.vcf_to_dataframe('test.vcf', fields='*', alt_number=2)
print(df2)

print(df2.query('DP > 10 and QUAL > 20'))


## preparation for larger vcf files --------------------------------------------
'''
callset = allel.read_vcf(vcf_path, fields=['numalt'], log=sys.stdout)

# largest number of alleles

numalt = callset['variants/numalt']
np.max(numalt)
'''

##----------------------------------
zarr_path = 'data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.zarr'
allel.vcf_to_zarr(vcf_path, zarr_path, group='22',
                  fields='*', alt_number=8, log=sys.stdout,
                  compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))

callset_h1k = zarr.open_group(zarr_path, mode='r')
callset_h1k
callset_h1k.tree(expand=True)
pos = allel.SortedIndex(callset_h1k['22/variants/POS'])

def plot_windowed_variant_density(pos, window_size, title=None):

    # setup windows
    bins = np.arange(0, pos.max(), window_size)

    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2

    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size

    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)

plot_windowed_variant_density(pos, window_size=100000, title='Variant density')
