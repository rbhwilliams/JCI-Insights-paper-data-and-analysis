# JCI-Insights-paper-data-and-analysis

This respository contains R analysis logs for RNA-Seq and related functional analysis from [Rousseau et al. 2016](https://insight.jci.org/articles/view/88178).

## Required data files and other downloads
Running this analysis requires a number of datafiles and R code or packages. See notes here for access:
1. The raw FASTQ files are available via NCBI BioProject Accession [PRJNA335539](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA335539). We have also provided a pre-processed version of the data [here] 
2. The ImmGen Consortium data is available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907. Please download the file 'GSE15907_RAW.tar' which contains the complete set of .CEL files from the study. These raw data are then processed using the R script in the file 'immgen.R' and results contained output in the file 'immgene-expression-data.RData'. This script is dependent on having installed several R/Bioconductor packages that are referenced in the first few line of 'immgen.R'. See https://www.bioconductor.org/install/ for more instructions.The file 'gene_assignment_from_immgen.txt' contains the gene memberships for 'coarse modules' defined by Jojic et al and is sourced from Supplementary Material of that paper.
3. Gene Ontology Enrichment Analysis. Builds on the R/Bioconductor package 'GO.db' with annotations sourced from the NCBI files 'gene2refseq' and 'gene2go', which are available from ftp://ftp.ncbi.nih.gov/gene/. See the README file on that page. The related R code is located in the file 'R.ontoTools.128.zip' and 'rw.ontoTools.R'. The code we have used here is *not recommended for general use* and was used here primarily as a vehicle for updating the now defunct Bioconductor package 'ontoTools'.
4. Miscellaneous: the R package 'identifier.mapping' available at https://github.com/rbhwilliams/identifier.mapping. The code used to compute the Betina statistic was obtained from: http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/EnrichmentProfilerDownload/EnrichmentScore.R

## Setting up and runnig the analysis
From the working directory, please make a directory called 'rfunctions' and place all .R files there. Unpack 'R.ontoTools.128.zip' and move it to rfunctions. Make a directory called 'datafiles' and place 'geneid_readcount_mm10', 'immgene-expression-data.RData' and 'gene_assignment_from_immgen.txt' there. Run R in the working directory and use the 'source' command on 'log.R' to run the entire analysis and produce all outputs as specified in the paper.
