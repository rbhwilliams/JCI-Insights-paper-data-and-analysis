# JCI-Insights-paper-data-and-analysis

This respository contains R analysis logs for RNA-Seq and related functional analysis from Rousseau et al. 2016. JCI Insights 1(15):e88178 https://insight.jci.org/articles/view/88178

The following files are either included in this respository or are available via the provided URLs:
1. The raw FASTQ files are available via https://submit.ncbi.nlm.nih.gov/subs/sra/SUB2597568/overview. A processed version of the data is made available in the file 'geneid_readcount_mm10.txt'.
2. The ImmGen Consortium data is available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907. Please download the file 'GSE15907_RAW.tar' which contains the complete set of .CEL files from the study. These raw data are then processed using the R script in the file 'immgen.R' and results contained output in the file 'immgene-expression-data.RData'. This script is depenedent on having installed seveal R/Bioconductor packages referenced in the first few line of 'immgen.R'. See https://www.bioconductor.org/install/ for more instructions.
3. The files 'gene2refseq' and 'gene2go' are available from ftp://ftp.ncbi.nih.gov/gene/. See the README file on that page.
The directory 'R' contains a number of files containing R code used in this study, mostly relating to enrichment analyses, and will be installed via source commands in each of the scripts. The code we have used here is not recommended for general use, as far more efficient implementations exist in the R/Bioconductor. The code was used here primarily as a vehicle for updating the now defunct Bioconductor package 'ontoTools'. 
4. The R package 'identifier.mapping' available at https://github.com/rbhwilliams/identifier.mapping
