#--updated analysis log for catheterisation paper--
#--started based on earlier logs by RW on 03 March 2016--

library(gdata)
library(DESeq2)
library(GO.db)
library(SparseM)
library(igraph)
library(graph)
library(limma)
library(sm)
library(AnnotationDbi)
library(mogene10sttranscriptcluster.db)
library(colorspace)
library(identifier.mapping)

#--SECTION: DIFFERENTIAL EXPRESSION ANALYSIS
#--read in the RefSeq annotated vesion of the expression data from HiSeq run 1
exp<-read.table(file="datafiles/geneid_readcount_mm10",header=T,sep="\t")
expCut<-as.matrix(exp[,2:7])
rownames(expCut)<-as.character(exp[,1])

#--there appears to be a number of rows with all-zero values...so remove them
allZeroInd<-which(rowSums(expCut)<1e-15)
expCut1<-expCut[-allZeroInd,]
#--and perform analysis using default approach and settings of DESeq2
colData<-data.frame(group=c(rep("naive",3),rep("catheter",3)))
rownames(colData)<-colnames(expCut1)
dds<-DESeqDataSetFromMatrix(countData=expCut1,colData=colData,design=~group)
dds$group<-relevel(dds$group,"naive")
dds1<-DESeq(dds)
res<-results(dds1)

#--reality check! everything should be on the line-of-identity..
plot(rowMeans(log2(expCut1[,4:6]))-rowMeans(log2(expCut1[,1:3])),res$log2FoldChange,las=1)
abline(0,1,lty=2)

#--we can use the summary(res,<adjpval>) to report summary results for d.e. analysis
#--see 'de-summary-results.txt' for data required for Table S1
#--examples...top 200 genes plotted in a MA-plot format
plot(res$baseMean,res$log2FoldChange,log="x",las=1)
points(res$baseMean[order(res$padj,decreasing=F)[1:200]],res$log2FoldChange[order(res$padj,decreasing=F)[1:200]],col=2,pch=16)

#--SECTION:ANNOTATIONS--
#--in shell...do the following...
system("grep ^10090 datafiles/gene2refseq > datafiles/gene2refseq10090")
system("cut -f1-2,4,16 datafiles/gene2refseq10090 > datafiles/gene2refseq10090-cut")
#--then load 'gene2refseq10090-cut' into R...
gene2refseq.10090.cut<-read.table(file="datafiles/gene2refseq10090-cut",header=F,sep="\t",stringsAsFactors=F)
gene2refseq.10090.cut<-gene2refseq.10090.cut[-which(gene2refseq.10090.cut[,1]!="10090"),]
gene2refseq.10090.cut<-gene2refseq.10090.cut[which(gene2refseq.10090.cut[,3]!="-"),]
#--the third column has .1 tags on the end of every refseq id...need to remove these...
gene2refseq.10090.cut.newCol3<-sapply(gene2refseq.10090.cut[,3],FUN=function(x){strsplit(x,".",fixed=T)[[1]][1]})
gene2refseq.10090.cut[,3]<-gene2refseq.10090.cut.newCol3
#--check and remove duplicated rows...'ur' denotes "/u/nique /r/ows"
gene2refseq.10090.cut.ur<-apply(gene2refseq.10090.cut[,2:3],1,FUN=paste,collapse="-")
gene2refseq.10090.cut.ur<-gdata::trim(gene2refseq.10090.cut.ur)
gene2refseq.10090.cut.ur.u<-unique(gene2refseq.10090.cut.ur)
gene2refseq.10090.cut.simple<-gene2refseq.10090.cut[match(gene2refseq.10090.cut.ur.u,gene2refseq.10090.cut.ur),2:4]

#--remove any genes that do not appear in our results...
gene2refseq.10090.cut.simple.exp<-gene2refseq.10090.cut.simple[which(gene2refseq.10090.cut.simple[,2]%in%rownames(res)==T),]
colnames(gene2refseq.10090.cut.simple.exp)<-c("geneid","refseq","symbol")

#--there can be one-one (mostly,hopefully), many-to-one and one-to-many and mant-many mappings
#--but looks like it's mostly going to be one-one...
#> apply(gene2refseq.10090.cut.simple.exp,2,FUN=function(x){length(unique(x))})
#geneid refseq symbol
#14226  14786  14226

#--now generate mapping between gene-ids and refseq-ids using my id2id mapping code...
#--we'll to cut the 3rd column of 'gene2refseq.10090.cut.simple.exp' and turn into a matrix...there is some complicated shit here with 'gdata::trim' (versions in IRanges and gdata)
gene2refseq.10090.cut.simple.exp.col12only<-cbind(gdata::trim(as.character(gene2refseq.10090.cut.simple.exp[,1])),gdata::trim(as.character(gene2refseq.10090.cut.simple.exp[,2])))
gene2refseq.10090.cut.simple.exp.idmap<-define.id2id.mapping.as.graph(gene2refseq.10090.cut.simple.exp.col12only)
gene2refseq.10090.cut.simple.exp.idmap.summ<-summarise.id2id.mapping(gene2refseq.10090.cut.simple.exp.idmap)

#--at this stage, we have to check that there is only one gene id that appears in each connected component
#--specifically, we can test this 'table(gene2refseq.10090.cut.simple.exp.idmap.summ$no.id1)'

#> table(gene2refseq.10090.cut.simple.exp.idmap.summ$no)
#2       3      4      5
#13701   494    27     4

#--good, so all gene ids have a single copy on this mapping, and only multiple refseq ids have multiple mappings!
#> table(gene2refseq.10090.cut.simple.exp.idmap.summ$no.id1)
#1
#14226

#> table(gene2refseq.10090.cut.simple.exp.idmap.summ$no.id2)
#1       2      3      4
#13701   494    27     4

#--for cases in which multiple refseq map to the same gene id, we can use the "choose the refseq with maximum average expression across all samples" filter
#--do so we'll need to calculate the mean expression...
res.meanexp.lg2.named<-log2(res$baseMean)
names(res.meanexp.lg2.named)<-rownames(res)

#--and then subselect one-to-many cases and filter (should probably turn this into a function, given I'm doing it twice, once for each replicate experiment)
multiCases<-gene2refseq.10090.cut.simple.exp.idmap$id2idCC[which(gene2refseq.10090.cut.simple.exp.idmap.summ$no>2)]
refseqExcludeFromMulti=NULL
new.gene2refseq.mapping=NULL
for(curcase in 1:length(multiCases))
{
 curgene<-intersect(V(multiCases[[curcase]])$name,gene2refseq.10090.cut.simple.exp.idmap$uid1)
 currefseq<-intersect(V(multiCases[[curcase]])$name,gene2refseq.10090.cut.simple.exp.idmap$uid2)
 currefseq.with.max.exp<-names(which.max(res.meanexp.lg2.named[currefseq]))[1]
 new.gene2refseq.mapping<-rbind(new.gene2refseq.mapping,c(geneid=curgene,refseq=currefseq.with.max.exp))
 refseqExcludeFromMulti<-c(refseqExcludeFromMulti,currefseq)
}

#--construct a set of one-one mappings from the one-to-one cases and filtered one-to-many cases just computed...
gene2refseq.10090.cut.simple.exp.single<-base::rbind(gene2refseq.10090.cut.simple.exp.col12only[!gene2refseq.10090.cut.simple.exp.col12only[,2]%in%refseqExcludeFromMulti,],new.gene2refseq.mapping)
gene2refseq.10090.cut.simple.exp.single.idmap<-define.id2id.mapping.as.graph(as.matrix(gene2refseq.10090.cut.simple.exp.single))
gene2refseq.10090.cut.simple.exp.single.idmap.summ<-summarise.id2id.mapping(gene2refseq.10090.cut.simple.exp.single.idmap)

#--next, we need to subset 'res' based on our final choices of included genes...
res.single<-res[gene2refseq.10090.cut.simple.exp.single[,2],]
#--and replace refseq identifiers with gene ids...
res.single.geneids<-res.single
rownames(res.single.geneids)<-gene2refseq.10090.cut.simple.exp.single[,1]

#--SECTION: GENE ONTOLOGY PREPARATION--
#--these functions are from Vincent Carey's ontoTools package (not current in R/Bioconductor).
#--with some modifications and extra code from me provided in 'rw.ontoTools.R'.
source(file="rfunctions/R.ontoTools.128/buildGOgraph.R")
source(file="rfunctions/R.ontoTools.128/namedSparse.R")
source(file="rfunctions/R.ontoTools.128/DAGtools.R")
source(file="rfunctions/R.ontoTools.128/graphExt.R")
source(file="rfunctions/R.ontoTools.128/makeSparseZero.R")
source(file="rfunctions/R.ontoTools.128/namedSparse.R")
source(file="rfunctions/R.ontoTools.128/ontoClasses.R")
source(file="rfunctions/R.ontoTools.128/ontoDepth.R")
source(file="rfunctions/rw.ontoTools.R")

#--calculate the ontology data..
goBPgraph<-buildGOgraph(useenv=GOBPPARENTS)
goBPDAG<-new("rootedDAG",root="GO:0008150",DAG=goBPgraph)
GOBP<-new("ontology",name="GOBF",version="this.one",rDAG=goBPDAG)
goBPamat<-clean.up.nS(accessMat(GOBP))
goBPTerms<-sapply(mget(nodes(goBPgraph),GOTERM),Term)

#--and import NCBI annotations...
system("grep ^10090 datafiles/gene2go > datafiles/gene2go.mus")
gene2go.mus.draft<-read.table(file="datafiles/gene2go.mus",header=F,sep="\t",as.is=T,quote="\"")
gene2go.mus.only<-gene2go.mus.draft[which(gene2go.mus.draft[,1]=="10090"),]
#--make a "bp" version...
gene2gobp.mus.no.nd<-remove.nd(filter.namespace(gene2go.mus.only,ns="Process"))
#--and match terms sets...
gene2gobp.mus.no.nd.matched<-filter.gene2go.for.goids(gene2gobp.mus.no.nd,GOBP)

#--and make the annotation matrices...
gene2gobp.mus.no.nd.matched.oomat<-make.oomat(gene2gobp.mus.no.nd.matched,GOBP)
gene2gobp.mus.no.nd.matched.covmat<-make.covmat(gene2gobp.mus.no.nd.matched.oomat,goBPamat)
exp1.gene2gobp.mus.no.nd.matched.covmat<-gene2gobp.mus.no.nd.matched.covmat
save(exp1.gene2gobp.mus.no.nd.matched.covmat,file="exp1.gene2gobp.mus.no.nd.matched.covmat.rda")

#--make an experiment specific annotation set...
gene2gobp.mus.no.nd.matched.covmat.v1<-cut.down.coverage.matrix.updated(gene2gobp.mus.no.nd.matched.covmat,rownames(res.single.geneids))
#--and calculate the inforamtion content of each term on this annotation set...
gene2gobp.mus.no.nd.matched.covmat.v1.termIC<-(-log2(getRowCount(t.nS(gene2gobp.mus.no.nd.matched.covmat.v1))/nrow(gene2gobp.mus.no.nd.matched.covmat.v1)))

#--SECTION: GO ENRICHMENT ANALYSIS--
#--we need to rank 'res.single.geneids' by (unadjusted p-value)
res.single.geneids.ranked<-res.single.geneids[order(res.single.geneids$pvalue,decreasing=F),]

#--ok, the names here are a bit silly I admit...'a' refers to analysis of up-regulated genes and 'b' refers to analysis of down-regulated genes..
a<-go.GO(gene2gobp.mus.no.nd.matched.covmat.v1,3,4,goBPTerms,res.single.geneids.ranked,c(100,200,500,1000,1500),"up")
a1<-tabulate.go.GO.output(a)
b<-go.GO(gene2gobp.mus.no.nd.matched.covmat.v1,3,4,goBPTerms,res.single.geneids.ranked,c(100,200,500,1000,1500),"down")
b1<-tabulate.go.GO.output(b)

#--SECTION: TABLE OUTPUTS FOR SUPPLMENTARY MATERIAL--
#--start with up-regulated genes...
a1.aug<-tabulate.go.GO.output.aug(a)
a1.aug1<-data.frame(GOid=rownames(a1.aug),a1.aug,term=goBPTerms[rownames(a1.aug)],stringsAsFactors=F)
a1.aug2<-a1.aug1[order(a1.aug1[,"adjP.100"],decreasing=T),]
write.table(a1.aug2,file="exp1-GOBP-enrichmentStatsSumm-up-aug.txt",sep="\t",col.names=T,row.names=F)

#--do same for down-regulated genes...
b1.aug<-tabulate.go.GO.output.aug(b)
b1.aug1<-data.frame(GOid=rownames(b1.aug),b1.aug,term=goBPTerms[rownames(b1.aug)],stringsAsFactors=F)
b1.aug2<-b1.aug1[order(b1.aug1[,"adjP.100"],decreasing=T),]
write.table(b1.aug2,file="exp1-GOBP-enrichmentStatsSumm-down-aug.txt",sep="\t",col.names=T,row.names=F)

#--note: the files 'exp-GOBP-enrichmentStatsSumm-up-aug.txt' and 'exp-GOBP-enrichmentStatsSumm-down-aug.txt' are imported into Excel, formatted and saved as "Table S2.xlsx"

#--SECTION: IMMGEN MODULE ANALYSIS--
source(file="rfunctions/immgen-functions.R")
#--I previously downloaded all .CEL files from GSE15907, normalised them using RMA and outputted everything into this .RData file (code not shown here)
load(file="datafiles/immgene-expression-data.RData")

#--we're using coarse modules from Immgen, the memberships are contained in the file 'gene_assignment_from_immgen.txt' which I sourced from the SOM of Jojic et al.
#--read this in, convert it to a list...
immgen<-read.table(file="datafiles/gene_assignment_from_immgen.txt",header=T,sep="\t")
immgen.gene2coarseList<-tapply(immgen$Gene,INDEX=immgen$Coarse.module,FUN=function(x){as.character(x)})
#--we need to convert these to gene ids...Immgen have done something really sloppy, which is they used  gene symbols with all-upper-case chars (for mouse!?!?!?), instead of official gene symbols
#--a fudge for this is to convert the (correct) gene symbols form the NCBI files I've sourced to upper case, before matching agains the coarse module gene memberships
immgen.gene2coarseList.geneid<-lapply(immgen.gene2coarseList,FUN=function(x,m){m[na.omit(match(x,toupper(m[,3]))),1]},m=gene2refseq.10090.cut.simple)
#--replace ids with expression values from 'res.single.geneids.ranked$log2FoldChange'
immgen.gene2coarseList.geneid.exp<-lapply(immgen.gene2coarseList.geneid,FUN=function(x,es){es$log2FoldChange[na.omit(match(x,rownames(es)))]},es=res.single.geneids.ranked)

#--do any of the modules have a mean log2 expression value that is stochastically different from zero?
#--so calculate one-sample t-test against mean of 0 (but we could also use a Mann-Whitney test here)
immgen.gene2coarseList.geneid.exp.notsmall<-immgen.gene2coarseList.geneid.exp[which(sapply(immgen.gene2coarseList.geneid.exp,length)>=3)]
immgen.gene2coarseList.geneid.exp.notsmall.osttest<-lapply(immgen.gene2coarseList.geneid.exp.notsmall,t.test)
immgen.gene2coarseList.geneid.exp.notsmall.osttest.res<-cbind(num=sapply(immgen.gene2coarseList.geneid.exp.notsmall,length),
                                                            estimate=sapply(immgen.gene2coarseList.geneid.exp.notsmall.osttest,FUN=function(x){x$estimate}),
                                                            rawp=-log10(p.adjust(sapply(immgen.gene2coarseList.geneid.exp.notsmall.osttest,FUN=function(x){x$p.value}),"bonf")))
rownames(immgen.gene2coarseList.geneid.exp.notsmall.osttest.res)<-names(immgen.gene2coarseList.geneid.exp.notsmall)
write.table(immgen.gene2coarseList.geneid.exp.notsmall.osttest.res,file="exp1-immgen.gene2coarseList.geneid.exp.notsmall.osttest.res.txt",sep="\t",col.names=T,row.names=T)

#--also need to generate data for density plots used in Figure S2--
immgen.gene2coarseList.geneid.exp.notsmall.smDensityDataList<-list()
for(s in (1:length(immgen.gene2coarseList.geneid.exp.notsmall)))
{
 curset<-names(immgen.gene2coarseList.geneid.exp.notsmall)[s]
 curres<-sm.density(immgen.gene2coarseList.geneid.exp.notsmall[[s]],display="none")
 immgen.gene2coarseList.geneid.exp.notsmall.smDensityDataList[[curset]]<-list(ep=curres$eval.points,estimate=curres$estimate)
}

#--for the part of the Immgen analysis that appears in Fig 1b and related SOM, we will need to compute the following sets of top-ranked genes...
#--up-regulated
top100.ps.up<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:100,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange>0)]),1])
top200.ps.up<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:200,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange>0)]),1])
top500.ps.up<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:500,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange>0)]),1])
#--down-regulated
top100.ps.down<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:100,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange<0)]),1])
top200.ps.down<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:200,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange<0)]),1])
top500.ps.down<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%intersect(rownames(res.single.geneids.ranked[1:500,]),rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange<0)]),1])

#--we need to construct the same lineage sets as in the website figure...
#--no easy way to do this, so I have copied entire strings from the Excel cells in Supplementary Table 2  of Jojic et al, then converted them to vectors using the function 'testfun'
lineageSetSampleTagList<-list(GN="GN.ARTH.BM,GN.ARTH.SYNF,GN.BL,GN.BM,GN.THIO.PC,GN.URAC.PC",
        MF="MF.103-11B+.LU,MF.103-11B+.SALM3.SI,MF.103-11B+.SI,MF.11CLOSER.SALM3.SI,MF.11CLOSER.SI,MF.169+11CHI.SLN,MF.BM,MF.II+480LO.PC,MF.II-480HI.PC,MF.LU,MF.MEDL.SLN,MF.MICROGLIA.CNS,MF.RP.SP,MF.SBCAPS.SLN,MF.THIO5.II+480INT.PC,MF.THIO5.II+480LO.PC,MF.THIO5.II-480HI.PC,MF.THIO5.II-480INT.PC",
        MO="MO.6C+II+.BL,MO.6C+II-.BL,MO.6C+II-.BM,MO.6C+II-.LN,MO.6C-II+.BL,MO.6C-II-.BL,MO.6C-II-.BM,MO.6C-IIINT.BL",
        DC="DC.103+11B+.SALM3.SI,DC.103+11B+.SI,DC.103+11B-.LU,DC.103+11B-.LULN,DC.103+11B-.LV,DC.103+11B-.POLYIC.LU,DC.103+11B-.SALM3.SI,DC.103+11B-.SI,DC.103-11B+.LULN,DC.103-11B+.LV,DC.103-11B+.POLYIC.LU,DC.103-11B+24+.LU,DC.103-11B+F4/80LO.KD,DC.4+.MLN,DC.4+.SLN,DC.4+.SP,DC.8+.MLN,DC.8+.SLN,DC.8+.SP,DC.8+.TH,DC.8-.TH,DC.8-4-11B+.MLN,DC.8-4-11B+.SLN,DC.8-4-11B+.SP,DC.8-4-11B-.MLN,DC.8-4-11B-.SLN,DC.8-4-11B-.SP,DC.IIHILANG+103+11BLO.SLN,DC.IIHILANG+103-11B+.SLN,DC.IIHILANG-103-11B+.SLN,DC.IIHILANG-103-11BLO.SLN,DC.LC.SK,DC.PDC.8+.MLN,DC.PDC.8+.SLN,DC.PDC.8+.SP,DC.PDC.8-.SP",
       B="B.FO.LN,B.FO.MLN,B.FO.PC,B.FO.SP,B.FRE.BM,B.FRE.FL,B.FRF.BM,B.GC.SP,B.MZ.SP,B.T1.SP,B.T2.SP,B.T3.SP,B1A.PC,B1A.SP,B1B.PC",
        NK="NK.49CI+.SP,NK.49CI-.SP,NK.49H+.SP,NK.49H-.SP,NK.B2M-.SP,NK.H+.MCMV1.SP,NK.H+.MCMV7.SP,NK.H-.MCMV1.SP,NK.MCMV1.SP,NK.MCMV7.SP,NK.SP",
        T4="T.4+8INT.TH,T.4.LN.BDC,T.4.PA.BDC,T.4.PLN.BDC,T.4FP3+25+.AA,T.4FP3+25+.LN,T.4FP3+25+.SP,T.4FP3-.SP,T.4MEM.LN,T.4MEM.SP,T.4MEM44H62L.LN,T.4MEM44H62L.SP,T.4NVE.LN,T.4NVE.MLN,T.4NVE.PP,T.4NVE.SP,T.4SP24-.TH,T.4SP24INT.TH,T.4SP69+.TH",
        T8="T.4INT8+.TH,T.8MEM.LN,T.8MEM.SP,T.8NVE.LN,T.8NVE.MLN,T.8NVE.PP,T.8NVE.SP,T.8SP24-.TH,T.8SP24INT.TH,T.8SP69+.TH",
        NKT="NKT.4+.LU,NKT.4+.LV,NKT.4+.SP,NKT.4-.LV,NKT.4-.SP,NKT.44+NK1.1+.TH,NKT.44+NK1.1-.TH,NKT.44-NK1.1-.TH",  GDT="TGD.SP,TGD.TH,TGD.VD1+24AHI.TH,TGD.VG1+VD6+.TH,TGD.VG1+VD6+.TH.ITKKO,TGD.VG1+VD6+.TH.TCRBKO,TGD.VG1+VD6+24AHI.TH,TGD.VG1+VD6+24AHI.TH.TCRBKO,TGD.VG1+VD6+24ALO.TH,TGD.VG1+VD6+24ALO.TH.TCRBKO,TGD.VG1+VD6-.TH,TGD.VG1+VD6-.TH.ITKKO,TGD.VG1+VD6-.TH.TCRBKO,TGD.VG1+VD6-24AHI.TH,TGD.VG1+VD6-24AHI.TH.ITKKO,TGD.VG1+VD6-24AHI.TH.TCRBKO,TGD.VG1+VD6-24ALO.TH,TGD.VG1+VD6-24ALO.TH.ITKKO,TGD.VG1+VD6-24ALO.TH.TCRBKO,TGD.VG2+.ACT.SP,TGD.VG2+.SP,TGD.VG2+.SP.TCRBKO,TGD.VG2+24AHI.E17.TH,TGD.VG2+24AHI.TH,TGD.VG2+24AHI.TH.ITKKO,TGD.VG2+24AHI.TH.TCRBKO,TGD.VG2+24ALO.TH,TGD.VG2+24ALO.TH.TCRBKO,TGD.VG2-.ACT.SP,TGD.VG2-.ACT.SP.TCRBKO,TGD.VG2-.SP,TGD.VG2-.SP.TCRBKO,TGD.VG2-24AHI.TH,TGD.VG2-24AHI.TH.TCRBKO,TGD.VG2-24ALO.TH.TCRBKO,TGD.VG3+24AHI.E17.TH,TGD.VG3+24ALO.E17.TH,TGD.VG5+.ACT.IEL,TGD.VG5+.IEL,TGD.VG5+24AHI.TH,TGD.VG5-.ACT.IEL,TGD.VG5-.IEL",
       SP="MLP.BM,MLP.FL,PROB.CLP.BM,PROB.CLP.FL,SC.CDP.BM,SC.CMP.BM,SC.GMP.BM,SC.LT34F.BM,SC.LTSL.BM,SC.LTSL.FL,SC.MDP.BM,SC.MEP.BM,SC.MPP34F.BM,SC.ST34F.BM,SC.STSL.BM,SC.STSL.FL",
        SC="EP.MECHI.TH,FI.MTS15+.TH,FI.SK,FRC.SLN,FRC.MLN,LEC.SLN,LEC.MLN,BEC.MLN,BEC.SLN,ST.31-38-44-.SLN")

#--reorder to match the Immgen figure website...
lineageSetSampleTagList1<-lineageSetSampleTagList[c("SP","B","DC","MF","MO","GN","T4","T8","NKT","GDT","SC","NK")]

#--convert and generate stuff we'll require for plots...
tf<-testfun(lineageSetSampleTagList1,immgenRmaMat)
tf<-tf[!duplicated(tf)]
#tc<-testcuts(tf)
#tcd<-diff(testcuts(testfun(lineageSetSampleTagList1,immgenRmaMat)))
#names(tcd)<-unique(names(tf))

#--we're using the analysis method of Benita (Ramnik Xavier's paper in Blood in 2013)
source(file="rfunctions/xavier-lab-method-code.R")
immgenRmaMat.v1<-immgenRmaMat[,tf]
colnames(immgenRmaMat.v1)<-names(tf)
immgenRmaMat.v1.enrichScore<-GenerateEnrichmentScore(immgenRmaMat.v1,names(tf))
immgenRmaMat.v1.enrichScore.summStat<-apply(immgenRmaMat.v1.enrichScore,2,FUN=quantile,probs=seq(0.01,0.99,0.01))

#--for our set of genes of interest, how many are above and below the 99th-%ile (and thus are consistent with infiltration)
#--generate random set of 100, 200 and 500 up and down regulated genes...
generate.random.probeset<-function(ps,N,B=100)
{
 res<-list()
 for(b in (1:B))
 {
  curpsset<-sample(ps,N,F)
  res[[b]]<-curpsset
 }
 return(res)
}

#--extract probesets from gene defined in our analysis
all.ps<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%rownames(res.single.geneids.ranked),1])
all.ps.up<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange>0)],1])
all.ps.down<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$log2FoldChange<0)],1])

#--and filter for those that are above minimum observed 'baseMean' for d.e. genes (so we get a realistic control set; see Fig S1 for details)
minMeanBaseForDeGenes<-min(res.single.geneids.ranked[which(res.single.geneids.ranked$padj<0.05),"baseMean"])
minMeanBaseGenes<-rownames(res.single.geneids.ranked)[which(res.single.geneids.ranked$baseMean>minMeanBaseForDeGenes)]
all.ps.minBaseMean<-unique(ps2EntrezList.no.NA.df[ps2EntrezList.no.NA.df[,2]%in%minMeanBaseGenes,1])

#--generate 100 random sets of genes sampling from 'all.ps.minBaseMean' of sizes 100,200 and 500...
#--we are computing separate sets for "up" and "down" to exactly match sample sizes in each sub-category...
top100.up.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top100.ps.up),100)
top200.up.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top200.ps.up),100)
top500.up.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top500.ps.up),100)
top100.down.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top100.ps.down),100)
top200.down.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top200.ps.down),100)
top500.down.randSet<-generate.random.probeset(all.ps.minBaseMean,length(top500.ps.down),100)

#--generate summary data for actual ("obs") and randomised ("rand") cases...
top100.up.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top100.ps.up)
top100.up.randSet.esCountMat<-sapply(top100.up.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)
top200.up.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top200.ps.up)
top200.up.randSet.esCountMat<-sapply(top200.up.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)
top500.up.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top500.ps.up)
top500.up.randSet.esCountMat<-sapply(top500.up.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)

#--and calculate summary statistics...
top100.down.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top100.ps.down)
top100.down.randSet.esCountMat<-sapply(top100.down.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)
top200.down.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top200.ps.down)
top200.down.randSet.esCountMat<-sapply(top200.down.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)
top500.down.obsSet.esCountVec<-gather.enrichment.score.statistics(immgenRmaMat.v1.enrichScore,top500.ps.down)
top500.down.randSet.esCountMat<-sapply(top500.down.randSet,FUN=gather.enrichment.score.statistics,esMat=immgenRmaMat.v1.enrichScore)

#--further outputs for paper...
#--make a table or Supplementary Data File from 'res.single.geneids.ranked'.
#--we'll need to add a couple of columns for Refseq id and Symbol...
rankedGeneSummTable<-data.frame(geneid=rownames(res.single.geneids.ranked),
                                    refseq=gene2refseq.10090.cut.simple.exp.single[match(rownames(res.single.geneids.ranked),gene2refseq.10090.cut.simple.exp.single[,1]),2],
                                    symbol=gene2refseq.10090.cut.simple[match(rownames(res.single.geneids.ranked),gene2refseq.10090.cut.simple[,1]),3],
                                    res.single.geneids.ranked,stringsAsFactors=F)
write.table(rankedGeneSummTable,file="exp1-rankedGeneSummTable.txt",sep="\t",row.names=F,col.names=T)

#--record versions and settings etc
sessionInfo()

#--and save workspace contents...
save.image(file="exp1.RData")

#--added on 14 July 2016--
#--plot for discussion with Molly et al...
m<-matrix(c(1,2,3,4,5,6),3,2)
layout(m,widths=c(1,1),heights=c(1,1,1))

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top100.up.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top100.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top100.up.randSet.esCountMat),rownames(top100.up.randSet.esCountMat))
lines(top100.up.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,11,c("top 100 de,logFC>0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"A",line=1.75)

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top200.up.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top200.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top200.up.randSet.esCountMat),rownames(top200.up.randSet.esCountMat))
lines(top200.up.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,25,c("top 200 de,logFC>0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"B",line=1.75)

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top500.up.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top500.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top500.up.randSet.esCountMat),rownames(top500.up.randSet.esCountMat))
lines(top500.up.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,60,c("top 500 de,logFC>0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"C",line=1.75)

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top100.down.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top100.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top100.down.randSet.esCountMat),rownames(top100.down.randSet.esCountMat))
lines(top100.down.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,11,c("top 100 de,logFC<0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"D",line=1.75)

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top200.down.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top200.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top200.down.randSet.esCountMat),rownames(top200.down.randSet.esCountMat))
lines(top200.down.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,25,c("top 200 de,logFC<0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"E",line=1.75)

par(mar=c(4.1,7.1,3.1,2.1))
matplot(t(apply(top500.down.randSet.esCountMat,1,fivenum)),type="l",ylim=c(0,max(top500.up.obsSet.esCountVec)),las=1,xaxt="n",lty=c(4,2,1,2,4),col="grey",lwd=c(1,1,2,1,1),ylab="Percentage of genes in top 1%\nof enrichment score distribution",xlab="Cell type")
axis(side=1,at=1:nrow(top500.down.randSet.esCountMat),rownames(top500.down.randSet.esCountMat))
lines(top500.down.obsSet.esCountVec,lwd=2,lty=1)
legend(6.5,60,c("top 500 de,logFC<0","random:max & min","random:1st & 3rd quartile","random:mean"),lty=c(1,4,2,1),lwd=c(2,1,1,2),col=c("black",rep("grey",3)),ncol=1,bty="n")
mtext(side=3,at=0.0,"F",line=1.75)

dev.print(pdf,file="exp1-Figure-S3-new.pdf",width=11,height=8)

