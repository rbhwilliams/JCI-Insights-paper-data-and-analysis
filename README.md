# JCI-Insights-paper-data-and-analysis

This respository contains R analysis logs for RNA-Seq and related functional analysis from [Rousseau et al. 2016](https://insight.jci.org/articles/view/88178).

## Required data files and other downloads
Running this analysis requires a number of datafiles and R code or packages. See notes here for access:
1. The raw FASTQ files are available via NCBI BioProject Accession [PRJNA335539](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA335539). We have also provided a pre-processed read count matrix with RefSeq ids indexed in rows and samples indexed in columns in the file [geneid_readcount_mm10](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/geneid_readcount_mm10) 
2. The ImmGen Consortium data is available from NCBI GEO Accession [GSE15097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907). Please download the file [GSE15907_RAW.tar](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15907&format=file) which contains the complete set of .CEL files from the study (**2.6Gb**). These raw data are then processed using the R script in the file [immgen.R](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/immgen.R) and results contained output in the file *immgene-expression-data.RData*. This script is dependent on having installed several R/Bioconductor packages that are referenced in the first few lines of [immgen.R](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/immgen.R). See https://www.bioconductor.org/install/ for more instructions.The file [gene_assignment_from_immgen.txt](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/gene_assignment_from_immgen.txt) contains the gene memberships for 'coarse modules' defined by [Jojic et al](http://www.nature.com/ni/journal/v14/n6/full/ni.2587.html) and was sourced from the Supplementary Material of that paper. Processing this data will take about one hour or so on a typical laptop. It is suggested to run *immgen.R* separately from the main analysis (use the *source* command in R to process *immgen.R*: see below for example) 
3. Gene Ontology Enrichment Analysis. This builds on the R/Bioconductor package [GO.db](http://bioconductor.org/packages/release/data/annotation/html/GO.db.html) with annotations sourced from the NCBI files *gene2refseq* and *gene2go*, which are available from via ftp://ftp.ncbi.nih.gov/gene/. See the README file on that page. The related R code is located in the file [R.ontoTools.128.zip](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/R.ontoTools.128.zip) and [rw.ontoTools.R](https://github.com/rbhwilliams/JCI-Insights-paper-data-and-analysis/blob/master/rw.ontoTools.R). The code we have used here is *not recommended for general use* and was used here primarily as a vehicle for updating the now defunct Bioconductor package *ontoTools*.
4. Other: you will also need the R package [identifier.mapping](https://github.com/rbhwilliams/identifier.mapping). The R [code](http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/EnrichmentProfilerDownload/EnrichmentScore.R) used to compute the Betina statistic was obtained from the Xavier lab website at MGH.  

## Setting up and running the analysis
1. From the working directory, make a directory called *rfunctions*.
2. Place all .R files in there. Unpack *R.ontoTools.128.zip* and also move them into *rfunctions*.
3. Also in the working directory, make a directory called *datafiles* and place *geneid_readcount_mm10*, *immgene-expression-data.RData* and *gene_assignment_from_immgen.txt* there.
4. Run R in the working directory and use:
```
> source(file="log.R")
```
to run the entire analysis and produce all outputs as specified in the paper.
