#--read in all CEL files and normalise using RMA--
#--requires install of R/Bioconductor packages oligo, AnnotationDbi and mogene10sttranscriptcluster.
#--run this script in the directory in which yuo uppack the .CEL files from the file 'GSE15907_RAW.tar'

library(oligo)
celFiles<-list.celfiles(full.names=F)
rawData<-read.celfiles(celFiles)
rmaData<-rma(rawData)

library(AnnotationDbi)
library("mogene10sttranscriptcluster.db")

#--we need versions that use both symbols and gene identifiers...
ps<-ls(mogene10sttranscriptclusterSYMBOL)
ps2SymbolList<-mget(ps,mogene10sttranscriptclusterSYMBOL)
ps2SymbolList.noNA<-ps2SymbolList[(!sapply(ps2SymbolList,is.na))]
ps2SymbolList.no.NA.df<-data.frame(ps=names(ps2SymbolList.noNA),symbols=unlist(ps2SymbolList.noNA),stringsAsFactors=F)

ps<-ls(mogene10sttranscriptclusterENTREZID)
ps2EntrezList<-mget(ps,mogene10sttranscriptclusterENTREZID)
ps2EntrezList.noNA<-ps2EntrezList[(!sapply(ps2EntrezList,is.na))]
ps2EntrezList.no.NA.df<-data.frame(ps=names(ps2EntrezList.noNA),symbols=unlist(ps2EntrezList.noNA),stringsAsFactors=F)

#--extract data..
immgenRmaMat<-assayData(rmaData)$exprs

#--save all this for use later..this file is used in the primary analysis
save(ps2SymbolList.no.NA.df,ps2EntrezList.no.NA.df,immgenRmaMat,file="immgene-expression-data.RData")
