#To extract genotypes for a given sample and a given gene, it's necessary to obtain the variant start
#and stop indices corresponding to the first and last variants within the gene.

## source: https://www.biostars.org/p/335077/
#---------------------------------------------------------------
# pick an arbitrary sample to work with
sample_idx = 42  # N.B., this is the 43rd sample, zero-based indexing
# load genotypes for the sample
gv = allel.GenotypeVector(callset[chrom]['calldata/GT'][:, sample_idx])

# setup some arrays to hold per-gene genotype counts for our sample of interest
import numpy as np
n_genes = len(genes_chr22)
n_hom_ref = np.zeros(n_genes, dtype=int)
n_het = np.zeros(n_genes, dtype=int)
n_hom_alt = np.zeros(n_genes, dtype=int)
n_variants = np.zeros(n_genes, dtype=int)

# iterate over genes
for i, (_, gene) in enumerate(genes_chr22.iterrows()):
    try:
        # locate data for this gene - this maps genomic coordinates onto variant start and stop indices
        loc_gene = pos.locate_range(gene.start, gene.end)
    except KeyError:
        # no data for the gene, leave counts as zero
        pass
    else:
        # extract genotypes for the gene
        gv_gene = gv[loc_gene]
        # compute genotype counts
        n_hom_ref[i] = gv_gene.count_hom_ref()
        n_het[i] = gv_gene.count_het()
        n_hom_alt[i] = gv_gene.count_hom_alt()
        # also store number of variants in the gene
        n_variants[i] = loc_gene.stop - loc_gene.start
