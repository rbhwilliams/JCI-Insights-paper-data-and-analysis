# JCI-Insights-paper-data-and-analysis

This respository contains R analysis logs for RNA-Seq and related functional analysis from Rousseau et al. 2016. JCI Insights 1(15):e88178 https://insight.jci.org/articles/view/88178

The following files are either included in this respository or are available via the provided URLs:
1. The raw FASTQ files are available via https://submit.ncbi.nlm.nih.gov/subs/sra/SUB2597568/overview. A processed version of the data is made available in the file 'geneid_readcount_mm10.txt'.
2. The ImmGen Consortium data is available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907. Please download the file 'GSE15907_RAW.tar' which contains the complete set of .CEL files from the study. These raw data are then processed using the R script in the file 'immgen.R' and results contained output in the file 'immgene-expression-data.RData'. This script is dependent on having installed several R/Bioconductor packages that are referenced in the first few line of 'immgen.R'. See https://www.bioconductor.org/install/ for more instructions.
3. Gene Ontology Enrichment Analysis. Builds on the R/Bioconductor package 'GO.db' with annotations sourced from the NCBI files 'gene2refseq' and 'gene2go', which are available from ftp://ftp.ncbi.nih.gov/gene/. See the README file on that page. The related R code is located in the directories HERE. The code we have used here is not recommended for general use, as far more efficient implementations exist in the R/Bioconductor, and was used here primarily as a vehicle for updating the now defunct Bioconductor package 'ontoTools'.
4. Miscellaneous: the R package 'identifier.mapping' available at https://github.com/rbhwilliams/identifier.mapping
5. The 'log.R' will execute the entire analysis and output figures and tables as produced in the paper. Pleae note you will have to set up directories for the above data objects as specified in the 'source' and 'load' commands in the first part of the file.
