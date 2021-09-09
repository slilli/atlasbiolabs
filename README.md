# atlasbiolabs

"""
@read_me:

    - hashing_ABO_permutations.py: dictionary for ABO permutations (A/B/AB/0) and genotypes (homozygous alternate/reference or heterozygous)
    - output.py: script to create .txt output, including samples, genotypes and bloodtypes
    - filter_relevant_snps.py: @here: 'rs8176719', 'rs8176746', 'rs8176747'
  	- .vcf must be in the same directory as this python script
	  - output file is created in the same directory as this python script
    
    !!!  we assume,that snps in input vcf files are on the forward strand, otherwise conversion is necessary !!!
"""

@required librarys:
 - https://pypi.org/project/scikit-allel/
 - numpy
 - pandas
 - sys

@sources:
  - https://www.biostars.org/p/42090/
  - https://www.snpedia.com/index.php/ABO_blood_group
  - http://alimanfoo.github.io/2017/06/14/read-vcf.html
  -  https://www.internationalgenome.org/data-portal/population/ASW
  -  https://www.biostars.org/p/302940/
  -  https://scikit-allel.readthedocs.io/en/stable/chunked.html#functions
  - https://blog.goldenhelix.com/alternate-allele-frequency-vcf-file-format/
  - https://pyvcf.readthedocs.io/en/v0.4.6/INTRO.html
  - https://scikit-allel.readthedocs.io/en/stable/genindex.html
  - https://stackoverflow.com/questions/5404665/accessing-elements-of-python-dictionary-by-index
