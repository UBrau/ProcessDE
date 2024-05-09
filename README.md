# ProcessDE
### Processing RNA-seq data for gene expression analysis

Starting from either isoform-level counts quantified with [_Salmon_](https://combine-lab.github.io/salmon/) (recommended), BAM files quantified with [_samtools_](http://www.htslib.org/doc/samtools-idxstats.html), or [_vast-tools_](https://github.com/vastgroup/vast-tools) expression counts, use the _exactTest_ functionality in the R package _edgeR_ to quantify each desired contrast. Replicates belonging to each sample type and contrasts are specified in CSV tables.
Genes are considered differentially expressed if (1) their expression, measured in counts per million (CPM), is at least a given threshold in at least one sample type of a contrast; (2) their log2-fold change is at least a given threshold, and (3) their FDR is lower than a given threshold. Defaults for thresholds can be changed optionally.

Note that if your experimental design is complex, this approach may not be appropriate. 

## Dependencies
R packages 'optparse', 'edgeR', 'gplots', and 'tximport'.

## Usage
`ProcessDE -h`

`ProcessDE -s <CSVFILE> -c <CSVFILE> --tx2gene <TSVFILE> --geneInfo <TSVFILE>`

## Input
- **Raw data files**: _Salmon_ transcript-level output files (recommended); output from _samtools idxstats_ quantification of BAM files; or _vast-tools_ expression output table
- **Sample table** in CSV format with columns _Sample_ and _Type_. Sample names must be unique. There must be replicates for at least one of the sample types in each contrast. In case of _Salmon_ counts, an additional _File_ column must be present that contains the path to the raw files.
- **Contrast table** in CSV format with columns _Experimental_ and _Control_. Each line defines a contrast, where the experimental and control types must match entries in the _Type_ column of the sample table.
- In case of _Salmon_ input: A tab-separated **file specifying transcript associations with genes** with, at least, columns _transcID_, _geneID_ and _geneName_, where transcID and geneID match the IDs used for _Salmon_ quantification.
- In case of _Salmon_ input: A tab-separated **gene information file** with, at least, columns _chrom_, _start_, _end_, _strand_, _geneID_, _geneName_, _biotype_.

## Output
- RPKM table per gene and sample
- Table with differentially expressed genes that fulfill the expression level, fold-change and FDR criteria; one per contrast
- Joint table with all differentially expressed genes from all contrasts
- MDS plot (a dimensionality reduction similar to PCA) of individual samples
- Correlation heatmap of individual sample pseudo-counts
- MA-plots (log-fold change vs. abundance) for each contrast, with differential genes highlighted
- Correlation heatmap of sample type RPKM
- Correlation heatmap of contrast log2-fold change (if more than one contrast)
- Scatter plots of log2-fold changes, comparing all contrasts (if more than one)
- Plot of the biological coefficient of variation from edgeR
- Tables of _tximport_ effective gene lengths and pseudocounts
- Optionally, folder with input files for GO analysis per contrast
- Log file

## Example input
Input files are provided in the _example_ folder based on _Salmon_ counts of RNA-seq samples from mouse neurons depolarized with KCl, Maze et al., _Neuron_ (2015), GEO: [GSE69807](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69807). The subfolder auxFiles contains a transcript to gene mapping file and a gene information file generated from the GENCODE vM21 annotation.

## Author
Please let me know if you have questions or encounter errors by emailing or raising an issue.

Ulrich Braunschweig, University of Toronto

[email](mailto:u.braunschweig@utoronto.ca)

## Reference
If using in published work, please cite the DOI of the release. The latest one is:
[![DOI](https://zenodo.org/badge/493020642.svg)](https://zenodo.org/badge/latestdoi/493020642)
