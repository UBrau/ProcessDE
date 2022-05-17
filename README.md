# ProcessDE
Processing RNAseq data for gene expression analysis

Starting from either isoform-level counts quantified with salmon preferred), BAM files, or vast-tools counts, use the exactTest functionality in the R package edgeR to quantify each desired contrast. Replicates belonging to each sample type and contrasts are specified in CSV tables.

## Dependencies
R packages 'optparse', 'edgeR', 'gplots', and 'tximport'.
