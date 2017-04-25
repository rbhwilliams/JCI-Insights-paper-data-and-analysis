#--functions for the paper "R software for clustering Gene Ontology annotations with application to human hepatitis transcriptomes"
#--some additional/modified functions added on 29 June 2015--

#---------------------------------------------------------------------------------
#--Functions for processing 'gene2go' files
#--The gene2go file can be sourced from the FTP site at NCBI/Entrez Gene database: 
#--ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#---------------------------------------------------------------------------------

#--two minor functions for processing 'gene2go'
#--a function to filter rows of 'gene2go' for the sub-ontolgy of interest
filter.namespace<-function(gene2go,ns=c("Component","Function","Process")){gene2go[which(gene2go[,8]==ns),]}
#--and a function to remove "ND" ("No biological Data available") annotations...
remove.nd<-function(gene2go){gene2go[which(gene2go[,4]!="ND"),]}

#--a function for checking consistency between GO tems in 'gene2go.gz' and BioC/GO.db
filter.gene2go.for.goids<-function(gene2go,ontology)
{
 #--'gene2go' is an R impornt from ncbi/gene2go.gz 
 #--'ontology' is an ontoTools object of class 'ontology'
 #--this function filter rows of 'gene2go' for terms that do not occur in 'ontology'
   
 res=NULL
 
 ind<-match(as.character(gene2go[,3]),nodes(DAG(rDAG(ontology))))
 res<-gene2go[which(is.na(ind)==F),] 
 
 return(res)

}

#---------------------------------------------------------------------------------
#--Functions for contructing object-ontology* matrices and coverage matrices 
#--*a binary gene-term matrix mapping primary annotations
#---------------------------------------------------------------------------------

#--a function to construct an object-ontology matrix from 'gene2go' and a corresponding ontology
make.oomat<-function(gene2go,ontology)
{
 #--a function to construct an 'OOmap' in ontoTools
 #--equivalent to the ooMap(x) component of an object of class object-ontology complex in ontoTools 
 #--'gene2go' is from ncbi/gene2go.gz 
 #--'ontology' is an ontoTools object of class 'ontology'
 #--'gene2go' should consistent with term set of ontology via application of 'filter.gene2go.for.goids'
 #--returns an BioC/ontoTools object of class "namedSparse"
 
 res=NULL
  
 #--define gene (gids) and terms (tids) in play...
 gids<-unique(as.character(gene2go[,2]))
 tids<-nodes(DAG(rDAG(ontology)))
 
 #--make a 'dummy' row to seed the sparse matrix construction--remove it after completion
 dud<-rep(0,length(tids))
 dud[1]<-1
 res<-as.matrix.csr(t(dud))
 #--use this to keep track of actual genes that are used in the matrix
 genes.to.use<-"dud"
 
 cat("Computing mappings for",length(gids),"gene ids...\n")
 for(i in (1:length(gids)))
 {
  if(i%%500==0)
  {
   cat(i)
  }
  #--extract anotations from gene2go for current gene id...
  tmp<-unique(as.character(gene2go[which(gene2go[,2]==gids[i]),3]))
  #--and test if there are any annotations, if so, store, otherwise get the next one
  if(length(tmp)>0)
  {
   tmpRow<-rep(0,length(tids))
   tmpRow[match(tmp,tids)]<-1
   res<-rbind(res,as.matrix.csr(t(tmpRow)))
   genes.to.use<-c(genes.to.use,gids[i])
  } 
 } 

 #--remove the "dummy" row
 res<-res[-1,]	    	    
 genes.to.use<-genes.to.use[-1]
 
 #--make into a namedSparse
 res<-new("namedSparse",mat=res,Dimnames=list(genes.to.use,tids))

 return(res)
 
}

#--a function to construct a coverage matrix from a object-ontology matrix and the accessibility matrix
make.covmat<-function(OOmap,accessMat)
{
 #--calculates a coverage matrix for an ontology and a set of primary annotations (contained in 'OOmap')
 #--coverage matrix is calculated as x+(x*y) were 'x'=object-ontology mapping matrix and 'y' is the accessibility-matrix.
 #--returns an BioC/ontoTools object of class "namedSparse"
 
 res=NULL

 a<-mat(OOmap)
 b<-a%*%mat(accessMat)
 res<-a+b
 res@ra<-pmin(1,res@ra)
 res<-new("namedSparse",mat=res,Dimnames=dimnames(OOmap))
 
 return(res)
 
}

#---------------------------------------------------------------------------------
#--Functions for working with namedSparse matrices 
#---------------------------------------------------------------------------------

#--a function to count the number of non-zero entries in a row of a namedSparse matrix, 'nS'
#--use 't(nS)' as an argument to count non-zero in columns of 'nS' 
getRowCount<-function(nS)
{
 #--counts non-zero rows elements of a namedSparse matrix
 #--'nS' is an BioC/ontoTools namedSparse
 
 res<-diff(mat(nS)@ia)
 names(res)<-dimnames(nS)[[1]]
 
 return(res)

} 


#--'cut.down.coverage.matrix':a function to subset rows of coverage matrix for genes of interest
cut.down.coverage.matrix<-function(nS,genes)
{
 #--'nS' is a namedSparse-matrix from ontoTools e.g. coverage matrix or object-ontology matrix. Should be genes-by-terms
 #--'genes' are a (unique) et of gene-ids that are a subset of dimnames(nS)[[1]]
 #--requires ontoTools,GO packages
 #--requires getRowCount function

 res=NULL
 ind=NULL
 colSum=NULL
 rowSum=NULL

 #--reduce to the geneset you're using...all of these should have at least one annootation, so we do not need to strip out any rows after this
 indRow<-na.omit(match(genes,dimnames(nS)[[1]]))
 res<-nS[indRow,]
 #--and compute colsum
 colSum<-getRowCount(t(res))
 #--and get rid of any terms that have no genes annotating to them...
 res<-res[,which(colSum>0)]  

 return(res)
 
}   

#--a function to remove columns from a coverage matrix to which all genes are annotated
remove.complete.columns<-function(covmat)
{
 #--a function to remove columns from a coverage matrix to which all genes are annotated
 #--'covmat' is a ontoTools coverage matrix
 
 res=NULL
 
 colSumInd<-which(getRowCount(t(covmat))==nrow(covmat))
 if(length(colSumInd)>0)
 {
  res<-covmat[,-colSumInd]
 }
 else
 {
  res<-covmat 
 }
  
 res
 
}

#--a function to remove genes that annotate to no terms
remove.non.annotated.rows<-function(covmat)
{
 #--a function that will remove any rows of a coverage that have zero rowSums
 #--'covmat' is a ontoTools coverage matrix
 
 res=NULL
 
 rowSumInd<-which(getRowCount(covmat)==0)
 if(length(rowSumInd)>0)
 {
  res<-covmat[-rowSumInd,]
 }
 else
 {
  res<-covmat
 }
 
 return(res)
	
}

#--'get.terms.for.a.gene': a function that will extract mapped terms from a object-ontology matrix (or coverage matrix) for a given gene...
get.terms.for.a.gene<-function(nS,geneid)
{
 #--will extract any terms (cols) in 'nS' that are annotated to a gene, 'geneid'
 #--'nS' should be a BioC/ontoTools namedSparse
 #--'geneid' is a character, specfiying a gene identifier (should be present in 'nS' i.e. an element of dimnames(nS)[[1]])
 #--returns a character vector of GO:id
 
 if(geneid%in%dimnames(nS)[[1]])
 {
  res<-dimnames(nS)[[2]][as.logical(as.matrix(mat(nS[geneid,])))]
  return(res)
 }
 else
 {
  cat("Gene not in OOmap\n")
 }
} 

#--'get.genes.for.a.term': a function that will extract genes that are annotated to a GO term
get.genes.for.a.term<-function(nS,goid)
{
 #--will extract any gene ids (rows) in 'nS' that are annotated to a term/GO:id, 'goid'
 #--'nS' should be a BioC/ontoTools namedSparse
 #--'goid' is a character, specfiying a GO term identifier (should be present in 'nS' i.e. an element of dimnames(nS)[[2]])
 #--returns a character vector of gene identifiers
  
if(goid%in%dimnames(nS)[[2]])
 {
  res<-dimnames(nS)[[1]][as.logical(as.matrix(mat(nS[,goid])))]
  return(res)
 }
 else
 {
  cat("GO:id not in OOmap\n")
 }
} 


#--'clean.up.nS': a function to remove 'zero' entries from the mat(nS)@ra slot of a spare matrix
#--the way ontoTools initialises sparse matrices leads to this problem
#--original version written by Mark Cowley while at UNSW (currenty: Garvan Institute of Medical Research, Sydney)
clean.up.nS<-function(nS)
{
 #--'nS' is a named sparse matrix from ontoTools

 #--local function
  rm.zeros.csr <- function(x)
  {
    if(length(x@ra) == 1 && x@ra[1] == 0)
        stop("Can't remove zero's since there's only one element")
    which <- rev( which(x@ra == 0) )
    if(length(which) == 0)
        return(x)
    x@ra <- x@ra[-which]
    x@ja <- x@ja[-which]
    ## update x@ia which is trickier
    for(i in which) {
        x@ia[x@ia > i] <- as.integer( x@ia[x@ia > i] - 1 )
    }
    return(x)
  }

  res=NULL

  res<-new("namedSparse",mat=rm.zeros.csr(mat(nS)),Dimnames=dimnames(nS))

  res

}

#-------------------------------------------------------------------------------------------------------------
#--SECTION 2: FUNCTIONS FOR CLUSTERING COVERAGE MATRICES
#-------------------------------------------------------------------------------------------------------------
#--note that this file has a modified version of 'cluster.covamt.extended'
cluster.coverage.matrix<-function(covmat,preproc=c("binary","cosine"),clustForm=c("single","complete"))
{
 #--function to cluster a coverage matrix
 #--uses hclust 
 #'covmat' is a ontoTools coverage matrix
 #'preproc' is the type of proprecessing you want to do:
 #'binary' is ubiqituous term removal with binary 
 #'cosine' is iaf-weighting and cosine similarity 

 #--requires functions--
 #--calc.perRow.cosines,getRowCount,iaf,calc.perRow.Norm2,remove.non.annotated.rows and remove.complete.columns
 
 #--performs preprocessing and 'hclust'-ering.
 #--has two outputs
 #--'res$clusters' is a list of vectors of gene-ids in each cluster
 #--'res$dendrogram' is a hclust object (useful for plots and further analysis, etc)
 #==in the local version, add the results of cutree directly...
 res=list(dendrogram=NULL,distMat=NULL)
   
 #--local-variables--'procCovMat' stores the processed coverage matrix
 procCovMat=NULL
 
 #--we'll need to process the coverage matrix for complete columns
 procCovMat<-remove.non.annotated.rows(remove.complete.columns(covmat))
 
 if(preproc=="binary")
 {
  #--make this dense, so we can add gene-ids as row.names
  tmp<-as.matrix(mat(procCovMat))
  rownames(tmp)<-dimnames(procCovMat)[[1]]
  colnames(tmp)<-dimnames(procCovMat)[[2]]
  res$distMat<-dist(tmp,"binary")
 }
 else if(preproc=="cosine") 
 {
  res$distMat<-as.dist((1-calc.perRow.cosines(iaf(procCovMat))))
 }
 
 #-generate dendrogram
 res$dendrogram<-hclust(res$distMat,clustForm)

 res

}

#--a function to form explicit clusters
make.cuts<-function(clusterObj,cutType,cut)
{
 #--a function to make cuts to a dendrogram and form explicit clusters
 #--arguments:
 #--'clusterObj' is an output of 'cluster.coverage.matrix'
 #--remaining args are for the function 'cutTreeLocal'
 #--'cutType' is the kind of cut you want to use in cuttree i.e. either "number" or "threshold"
 #--'cut' is a numeric scalar should correspond to 'cutType' 
 #--outputs:
 #--'clusters' is a list, whose elements contain the members of each group.
 #--'cutRes' stores the output of 'cutTreeLocal'
 #--"cutType' and 'cut' pass input args of same name
 
 res<-list(clusters=list(),cutRes=NULL,cutType=cutType,cut=cut)
 
 #--local variables--'ctr' stores the cut tree results; 'groups' stores the cluster-id using left-right dendrogram ordering
 ctr=NULL
 groups=NULL
 
 res$cutRes<-cutTreeLocal(hc=clusterObj$dendrogram,cutType,cut)
 #--and extract cluster using the dendrogram-ordering of groups...
 groups=unique(res$cutRes$Groups)
 #--and output list...
 for(i in (1:length(groups)))
 {
  res$clusters[[i]]<-as.vector(res$cutRes$Labels[which(res$cutRes$Groups==groups[i])])
 } 
 
 res
 
}   

calc.perRow.Norm2<-function(mat)
{
 #--fucntion to calculate the 2-norm (Euclidean) of each row in a matrix ('mat')
 #--returns vector of norms

 res<-apply(mat,1,FUN=function(x){sqrt(x%*%x)})

res

}

calc.perRow.cosines<-function(mat)
{
 #--function to calculate cosines between pairs of row vectors in a matrix ('mat')
 #--done the dumbo-way, full-symmetric matrix computed. Convert to dist afterwards if need be...
 
 res=NULL
 #pre-calc 2-norm for all rows
 perRow.2norm<-calc.perRow.Norm2(mat)

 for(i in (1:nrow(mat)))
 {
  res<-cbind(res,apply(mat,1,'%*%',mat[i,])/(perRow.2norm*perRow.2norm[i]))
 }

res
  
}

iaf<-function(covmat)
{
 #--function to weight a coverage matrix using inverse annotation frequency scheme
 #--note that ubiqutious terms have been removed
 #--'covmat' is an ontoTools (namedSparse) coverage matrix
 #--uses function 'getRowCount'

 res=NULL
 weights=NULL
 colSum=NULL
 colSS=NULL
 matlabels=NULL
 
 #--copy dimnames...
 matlabels<-dimnames(covmat)
 
 #compute weights and substitute with an element wise multiplication
 weights<-log(nrow(covmat)/getRowCount(t(covmat)))
 covmat<-t(as.matrix(mat(t(covmat)))*weights)
 
 #calculate normalisation factors for each column
 #we need to sum the squares of the each column element
 colSS<-sqrt(apply(covmat^2,2,sum))
 #and divide covmat column-wise
 res<-t(t(covmat)/colSS)
 rownames(res)<-matlabels[[1]]
 colnames(res)<-matlabels[[2]]
 
 res

}

cutTreeLocal<-function(hc,type,c)
{ 
 #--function to return results of a cutree operation
 #--pass usual arguments to cutree e.g. h or k
 #--hc is an hclust object
 #--labels and defined-clusters are returned according to hc$order.
 res=NULL 
 hc.groups=NULL

 if(type=="number")
 {
  hc.groups<-cutree(tree=hc,k=c,h=NULL)
 }
 else if(type=="threshold")
 {
  hc.groups<-cutree(tree=hc,k=NULL,h=c)
 }
 
 res<-data.frame(hc$labels[hc$order],hc.groups[hc$order],row.names=NULL)
 names(res)<-c("Labels","Groups") 
 
 res
 
}

#-------------------------------------------------------------------------------------------------------------
#--SECTION 2: FUNCTIONS FOR INTERPRETATING COVERAGE-MATRIX CLUSTERS
#-------------------------------------------------------------------------------------------------------------
#--this function is used to collect basic information on the characteristics of a cluster
#--it can be used in specialised plotting functions.
tag.clusters.by.terms<-function(covmat,clustRes,ds,scale)
{
 #---a function to extract the minimum term depth to which all genes map
 #--'covmat' is an ontoTools coverage matrix which has been clustered
 #--'clustRes' is an output from make.cuts
 #--'ds' is a "depth structure" object from ontoTools i.e. depthStruct(goBP170DAG)
 #--'scale' is a factor between 0 and 1. scale*(#(genes)) is a threshold used to select terms that characterise the cluster
 #--values for in the res$summary  matrix
 #--'d.comm' is the minimum depth term for which all genes map to is selected
 #--'d.min' is the lowest term depth in the cluster (should be zero)
 #--'d.max' is the highest term depth in the cluster
 #--'N.terms' is the number of terms in the cluster
 #--'N.genes' is the number of genes in the cluster
 #--also returns the actually go-ids that map to all genes in a cluster in the list res$termsList
 
 res=list(termsList=list(),summary=NULL)
 tmp=NULL
 terms=NULL
 
 res$summary<-array(NA,c(length(clustRes$clusters),5))
 
 for(i in (1:length(clustRes$clusters)))
 {
  #--extract the coverage matrix for the genes in each cluster and sum across rows
  tmp<-cut.down.coverage.matrix(covmat,clustRes$clusters[[i]])
  terms<-dimnames(tmp)[[2]][which(getRowCount(t(tmp))>=(scale*length(clustRes$clusters[[i]])))]
  res$summary[i,]<-c(max(unlist(lapply(terms,ds$tag2depth))),
						min(unlist(lapply(dimnames(tmp)[[2]],ds$tag2depth))),
						max(unlist(lapply(dimnames(tmp)[[2]],ds$tag2depth))),
						length(dimnames(tmp)[[2]]),length(clustRes$clusters[[i]]))
  res$termsList[[i]]<-terms[which(lapply(terms,ds$tag2depth)==max(unlist(lapply(terms,ds$tag2depth))))]
 }
  
 colnames(res$summary)<-c("d.common","d.min","d.max","N.terms","N.genes")  
  
 res
 
}

#--and a function to compute the information content of all annotated terms
calc.term.ic<-function(covMat,root)
{
 tmp<-getRowCount(t(covMat))
 
 ic<-(-log2(tmp/tmp[which(dimnames(covMat)[[2]]==root)]))

 ic
 
}

#-------------------------------------------------------------------------------------------------------------
#--SECTION 3: FUNCTIONS FOR VISUALISING CLUSTERING OUTPUT	
#-------------------------------------------------------------------------------------------------------------
extract.terms<-function(goidList,terms)
{
 #--'goidList' is a list, whose elements contain one or more GO-ids
 #--'terms' is a named vector, whose elements are term names, and whose names are GOids
 
 res=NULL
 tmp=NULL

 for(i in (1:length(goidList)))
 {
  tmp=NULL
  tmp<-as.character(terms[match(goidList[[i]],names(terms))])
  res<-c(res,paste(tmp,collapse="+"))
 }
 
 return(res)
 
}

#---now we need to define the points along the dendrogram where we will place the terms and 'level-points'
extract.place.points.on.dendro<-function(cutreeRes)
{
 #--'cutreeRes' is the output from cutree
 #--samples in 'cutreeRes' are listed as per left-right ordering on dendrogram
 #--find the unique groups (these will be in the same order as on the dendrogram
 #--'res' is the x-ords of the points along the 1:N axis for the N samples 
 res=NULL
 ind=NULL

 ind<-1:nrow(cutreeRes)
  
 groups<-unique(cutreeRes[,2])
  
 for(i in (1:length(groups)))
 {
  res<-c(res,mean(ind[which(cutreeRes[,2]==groups[i])]))
 }

 res

}

plot.clustering.results<-function(clustRes,cutRes,covmat,depth.structure,geneNames,terms)
{
 #--'clustRes' is an output from the function 'cluster.coverage.matrix'
 #--'cutRes' is an output from the function 'make.cuts'
 #--'covmat' is the ontoTools coveragte matrix that was used for the clustering
 #--'depth.structure' is e.g. depthStruct(goBP170DAG)
 #--'geneNames is a vector alternative node names for the dendrogam e.g.HUGO, etc.
 #--'terms' is a vector of named GO terms; names should be corresponding GOids.

 nf<-layout(matrix(c(1,2),2,1,byrow=TRUE),c(lcm(35),lcm(35)),c(lcm(6),lcm(17)))
  
 tagRes=NULL

 #--first, plot the dendrogram...
 par(mar=c(0.5,4,0.5,1))
 plot(clustRes$dendrogram,main="",axes=F,ann=F,labels=geneNames,cex=0.7)
 #--and overlay the clusters
 if(cutRes$cutType=="threshold")
 {
  my.rect.hclust(clustRes$dendrogram,k=NULL,h=cutRes$cut,border=1)
 }
 else
 {
  my.rect.hclust(clustRes$dendrogram,h=NULL,k=cutRes$cut,border=1)
 }

 #--compute the tags...
 tagRes<-tag.clusters.by.terms(covmat,cutRes,depth.structure,scale=1)
 #--now, compute the mid-points for the names, and the names themselves
 tagTerms<-extract.terms(tagRes[[1]],terms)
 tagPoints<-extract.place.points.on.dendro(cutRes$cutRes)
 axis(side=1,at=tagPoints,labels=tagTerms,las=2,cex.axis=0.7)

}

#--a modified version of rect.hclust
my.rect.hclust<-function(tree,k=NULL,which=NULL,x=NULL,h=NULL,border=2,cluster=NULL) 
{
    if (length(h) > 1 | length(k) > 1) 
        stop("'k' and 'h' must be a scalar")
    if (!is.null(h)) {
        if (!is.null(k)) 
            stop("specify exactly one of 'k' and 'h'")
        k <- min(which(rev(tree$height) < h))
        k <- max(k, 2)
    }
    else if (is.null(k)) 
        stop("specify exactly one of 'k' and 'h'")
    if (k < 2 | k > length(tree$height)) 
        stop(gettextf("k must be between 2 and %d", length(tree$height)), 
            domain = NA)
    if (is.null(cluster)) 
        cluster <- cutree(tree, k = k)
    clustab <- table(cluster)[unique(cluster[tree$order])]
    m <- c(0, cumsum(clustab))
    if (!is.null(x)) {
        if (!is.null(which)) 
            stop("specify exactly one of 'which' and 'x'")
        which <- x
        for (n in 1:length(x)) which[n] <- max(which(m < x[n]))
    }
    else if (is.null(which)) 
        which <- 1:k
    if (any(which > k)) 
        stop(gettextf("all elements of 'which' must be between 1 and %d", 
            k), domain = NA)
    border <- rep(border, length.out = length(which))
    retval <- list()
    for (n in 1:length(which)) {
        rect(m[which[n]] + 0.66, par("usr")[3], m[which[n] + 
            1] + 0.33, mean(rev(tree$height)[(k - 1):k]), border = border[n],lty=2)
        retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
    }
    invisible(retval)
}


augmentDot<-function(clusterGenes,refGenes,covmat,oomap,goGraph,goTerms)
{
 #--function to generate an augmented dot file for rendering in GraphViz
 #--'clusterGenes' is a vector of gene.ids e.g. from 'makeCuts'
 #--'refGenes' set of gene symbol-named gene.ids--NOTE that clusterGenes MUST be a subset of refGenes. 'refGenes' should have HUGO symbols as 'names'
 #--'covmat' is an ontoTools coverage matrix
 #--'oomap' is an ontoTools OOmap matrix
 #--'goGraph' is a graph of the ontology you are using
 #--'goTerms' is a vector GO terms, named by their corresponding GOids
	
 require(ontoTools)
	
 #requires function: extract.GO.for.geneList
	
 #declarations
 G=NULL
 G.dot=NULL
 geneNames=NULL
 genesNodesString=NULL
 genesEdgesString=NULL
 terms=NULL
 res=NULL

 #--subset all terms that are fully annotated to genes in a 'clusterGenes' (and format as .dot)	 
 G<-subGraph(dimnames(cut.down.coverage.matrix(covmat,clusterGenes))[[2]],goGraph)
 G.dot<-toDot(G)
	
 #--first, substitute solid lines for all GO connections
 G.dot<-gsub("dir=back",replacement="dir=back,style=solid",G.dot) 
	
 #--get HUGO-symbols for labels...
 genesNames<-names(refGenes)[match(clusterGenes,refGenes)]
	
 #--construct string for all gene-name nodes...
 for(i in (1:length(clusterGenes)))
 {
  genesNodesString<-paste(genesNodesString,paste("\"",as.character(genesNames[i]),"\"","[shape=box]",";\n", sep=""),collapse=" ")
 }
	
 #--generate edges...will need to access data in OOmap(OOC)..
 for(i in (1:length(clusterGenes)))
 {
  terms<-get.terms.for.a.gene(oomap,clusterGenes[i])
  cat(genesNames[i],"\t",terms,"\n")
  for(j in (1:length(terms)))
  {
   cat(paste("edge [dir=back,style=dashed] ","\"",terms[j],"\"","->","\"",genesNames[i],"\"",";\n",sep=""),"\n")
   genesEdgesString<-paste(genesEdgesString,paste("edge [dir=back,style=dashed] ","\"",terms[j],"\"","->","\"",genesNames[i],"\"",";\n",sep=""),sep=" ")
  } 
 }
	
 #--split up original string by occurence of 'edge'
 G.dot.split<-strsplit(G.dot,"edge")
 #--add the 'nodes' bit to fiveG.split[[1]][1], followed by the new edges.../
 G.dot.new<-paste(G.dot.split[[1]][1],genesNodesString,genesEdgesString,sep=" ")
 #--loop through all remaining elements and paste back in--make sure to add 'edge' to each...
 for(i in (2:length(G.dot.split[[1]])))
 {
  G.dot.new<-paste(G.dot.new,paste("edge",G.dot.split[[1]][i],sep=" "))
 }
	
 #--need a function that substitutes GO term names instead of ids in a dot file--script it here..
 for(i in (1:length(nodes(G))))
 {
  G.dot.new<-gsub(nodes(G)[i],goTerms[nodes(G)[i]],G.dot.new)
 }
	
 #--output string...
 return(G.dot.new)	
	
}

#--UPDATED FUNCTIONS ON 29 JUNE 2015--
#--and we might need one that transposes a nS...
t.nS<-function(nS)
{
 m<-nS@mat
 mt<-SparseM::t(m)
 res<-new("namedSparse",mat=mt,Dimnames=list(dimnames(nS)[[2]],dimnames(nS)[[1]]))
 return(res)
}

cut.down.coverage.matrix.updated<-function(nS,genes)
{
 #--'nS' is a namedSparse-matrix from ontoTools e.g. coverage matrix or object-ontology matrix. Should be genes-by-terms
 #--'genes' are a (unique) et of gene-ids that are a subset of dimnames(nS)[[1]]
 #--requires ontoTools,GO packages
 #--requires getRowCount function
    
 res=NULL
 #--reduce to the geneset you're using...all of these should have at least one annotation, so we do not need to strip out any rows after this
 indRow<-na.omit(match(genes,dimnames(nS)[[1]]))
 res<-nS[indRow,]
 #--and compute colsum
 colSum<-getRowCount(t.nS(res))
 #--and get rid of any terms that have no genes annotating to them...
 res<-res[,which(colSum>0)]
 return(res)
}

filter.coverage.matrix<-function(covmat,icData,ic.min=2,ic.max=3)
{
 #--returns a coverage matrix filtered against term information content range
 #--from the coverage matrix, remove terms that are outside the information content window
 termsToUse<-names(icData)[intersect(which(icData>=ic.min),which(icData<ic.max))]
    
 #--and then subset this coverage matrix...
 ccovmat<-cut.down.coverage.matrix.updated(covmat[,termsToUse],dimnames(covmat)[[1]])
    
 return(ccovmat)
}

#--convert to a full matrix...
convert2full<-function(covmatCsr)
{
 res<-as.matrix(covmatCsr@mat)
 rownames(res)<-dimnames(covmatCsr)[[1]]
 colnames(res)<-dimnames(covmatCsr)[[2]]
 return(res)
}

extract.nonzero.elements<-function(memmat,key)
{
 #--assume we're extracting from rows, transpose memmat if you want to extract from columns...
 res<-colnames(memmat)[which(memmat[key,]>0)]
 return(res)
}

convert.gene2pwMat.to.pw2geneList<-function(gene2pwMat)
{
 res<-list()
 for(pw in (1:ncol(gene2pwMat)))
 {
  curpw<-colnames(gene2pwMat)[pw]
  res[[curpw]]<-extract.nonzero.elements(t(gene2pwMat),curpw)
 }
 return(res)
}

process.pathway.list.against.data<-function(pwl,genes,min=5)
{
 #--'pwl' is an output from 'process.pathways.to.list'
 #--'genes' is a set of genes that are expressed/etc (symbols)
 #--'min' is an integer, specifying the minimum number of genes that can be in a pathway.
 res=NULL
 res<-lapply(pwl,FUN=intersect,y=genes)
 len<-sapply(res,length)
 res<-res[which(len>=min)]
 return(res)
}

fishers.test.for.a.pathway<-function(tst,bkg,pathwayMemberGenes)
{
 #--'tst' and 'bkg' are vectors of genes, and 'pathwayMemberGenes' contains genes in each pathway
 #--make a contingency table and return it
    
 #--categorise...
 A<-length(intersect(tst,pathwayMemberGenes))
 B<-length(setdiff(tst,pathwayMemberGenes))
 C<-length(intersect(bkg,pathwayMemberGenes))
 D<-length(setdiff(bkg,pathwayMemberGenes))
    
 #--make matrix...
 res<-matrix(c(A,B,C,D),2,2)
 rownames(res)<-c("in","out")
 colnames(res)<-c("de","nde")
 return(res)
}

run.pw.analysis<-function(pwl,tst,bkg)
{
 res<-list(sigContab=NULL,rawp=NULL,adjp=NULL,cellASets=NULL)
 contabList<-lapply(pwl,FUN=fishers.test.for.a.pathway,tst=tst,bkg=bkg)
 countCellA<-sapply(contabList,FUN=function(x){x[1,1]})
 res$rawp<-sapply(contabList,FUN=function(x){fisher.test(x,alternative="greater")$p.value})
 adjp<-sort(p.adjust(res$rawp,"BH"),decreasing=F)
 res$adjp<-adjp[which(countCellA[names(adjp)]>0)]
 res$sigContab<-contabList[names(res$adjp)]
 res$cellASets<-lapply(pwl[names(res$adjp)],FUN=intersect,y=tst)
 return(res)
}

compute.DAVID.output<-function(pwl,tst,universe,terms)
{
 #--computes GO-enrichment statistics a la DAVID :)
 res=NULL
 bkg<-setdiff(universe,tst)
 res<-list()
 res$contabList<-lapply(pwl,FUN=fishers.test.for.a.pathway,tst=tst,bkg=bkg)
 res$numAnnoInTst<-sapply(res$contabList,FUN=function(x){x["in","de"]})
 res$propAnnoInTst<-sapply(res$contabList,FUN=function(x){x["in","de"]/sum(x[,"de"])})
 res$propDeInAnno<-sapply(res$contabList,FUN=function(x){x["in","de"]/sum(x["in",])})
 res$expectedAnnoProp<-sapply(res$contabList,FUN=function(x){sum(x["in",])/sum(x)})
 res$obsPropVsexpProp<-res$propAnnoInTst/res$expectedAnnoProp
 res$pval<-sapply(res$contabList,FUN=function(x){fisher.test(x,alternative="greater")$p.value})
 res$terms<-terms[names(res$pval)]
 return(res)
}

make.table<-function(enrichmentList)
{
 res=NULL
 res<-data.frame(numAnnoInTst=enrichmentList$numAnnoInTst,
 propAnnoInTst=enrichmentList$propAnnoInTst,
 propDeInAnno=enrichmentList$propDeInAnno,
 expectedAnnoProp=enrichmentList$expectedAnnoProp,
 obsPropVsexpProp=enrichmentList$obsPropVsexpProp,
 rawp=enrichmentList$pval,
 adjp=p.adjust(enrichmentList$pval,"BH"),
 terms=enrichmentList$terms)
 res<-res[order(res$rawp,decreasing=F),]
 return(res)
}

run.DAVID.using.CAT<-function(rankedResTable,pwl,universe,term,catRun=c(100,200,300,500,1000),type=c("up","down"))
{
 #--take each of the the top 1:N genes, with the N-series defined in the vector 'catRun'
 #--run a DAVID analysis, return each result a list element...
 res<-list()
    
 #--run through each set of genes...
 for(n in (1:length(catRun)))
 {
  #--we need to define the local 'tst' argument...record 'curn' as well due to drop out
  if(type=="up")
  {
   curgenes<-intersect(rownames(rankedResTable)[1:catRun[n]],rownames(rankedResTable)[which(rankedResTable$log2FoldChange>0)])
  }
  else
  if(type=="down")
  {
   curgenes<-intersect(rownames(rankedResTable)[1:catRun[n]],rownames(rankedResTable)[which(rankedResTable$log2FoldChange<0)])
  }
  curtst<-intersect(curgenes,universe)
  curn<-length(curtst)
  print(curn)
  res[[n]]<-compute.DAVID.output(pwl,curtst,universe,term)
 }
 return(res)
}

make.term.wise.distance.matrix<-function(covmat)
{
 #--convert to (full) matrix format and cluster columns...
 covmatfull<-SparseM::as.matrix(covmat@mat)
 rownames(covmatfull)<-dimnames(covmat)[[1]]
 colnames(covmatfull)<-dimnames(covmat)[[2]]
 res<-as.matrix(dist(t(covmatfull),"binary"))
 return(res)
}

scale.Dmat<-function(Dmat)
{
 #--'Dmat' is a symmetric (distance) matrix of class 'matrix'
 res<-(-0.5*Dmat^2)
 return(res)
}

double.center.data.matrix<-function(m)
{
 rM<-rowMeans(m)
 cM<-colMeans(m)
 mM<-mean(as.vector(m))
    
 res<-matrix(NA,nrow(m),ncol(m))
 for(i in (1:nrow(m)))
 {
  for(j in (1:ncol(m)))
  {
   res[i,j]<-m[i,j]-rM[i]-cM[j]+mM
  }
 }
    
 rownames(res)<-rownames(m)
 colnames(res)<-colnames(m)
    
 return(res)
    
}

#--convert cartesian to polar...
c2p<-function(mat)
{
 #--returns angle to x-axis in degrees, along with radius of vector
 r<-sqrt(mat[,1]^2+mat[,2]^2)
 phi<-atan2(mat[,2],mat[,1])*(180/pi)
 ind<-which(phi<0)
 phid<-phi
 phid[ind]<-(phi[ind]+360)
 res<-cbind(phid,r)
 rownames(res)<-rownames(mat)
 return(res)
}

#--write a function to process everything...
go.GO<-function(covmat,icMin,icMax,goterms,rankedGeneTable,catRun,type=c("up","down"))
{
 #--do any entire GO enrichment analysis and visualisation!
 res<-list(icMin=icMin,icMax=icMax,termIC=NULL,covmat=NULL,catRun=catRun,enrichmentRes=NULL)
 
 #--calculate term information content...
 print(paste("Calculating information content of all",nrow(covmat),"terms",sep=" "))
 res$termIC<-(-log2(getRowCount(t.nS(covmat))/nrow(covmat)))
 
 #--process coverage matrix for IC thresholds
 res$covmat<-filter.coverage.matrix(covmat,res$termIC,icMin,icMax)
 print(paste("Coverage matrix extracted for",ncol(res$covmat),"terms with information content between",icMin,"and",icMax,sep=" "))

 #--generate lists of terms...
 print("Generate data for enrichment analysis...")
 covmat.full<-convert2full(res$covmat)
 covmat.full2list<-convert.gene2pwMat.to.pw2geneList(covmat.full)
 #--and compute DAVID results...
 print("...and compute enrichment results")
 res$enrichmentRes<-run.DAVID.using.CAT(rankedGeneTable,covmat.full2list,dimnames(res$covmat)[[1]],goterms,catRun,type)

 return(res)
 
}

tabulate.go.GO.output<-function(goGO)
{
 #--take an output from 'go.GO' and make a table
 goGO.adjpSumm<-sapply(goGO$enrichmentRes,FUN=function(x){-log10(p.adjust(x$pval))})
 goGO.adjpSummTable<-data.frame(goGO.adjpSumm,goBPTerms[rownames(goGO.adjpSumm)])
 colnames(goGO.adjpSummTable)<-c(as.character(goGO$catRun),"Term")
 goGO.adjpSummTable1<-goGO.adjpSummTable[order(rowMeans(goGO.adjpSumm),decreasing=T),]
 return(goGO.adjpSummTable1)
}

pcao.on.goGO<-function(goGO)
{
 #--take a IC-subsetted coverage matrix and analyse using PCoA
 dMat<-make.term.wise.distance.matrix(goGO$covmat)
 dMat.scale<-scale.Dmat(dMat)
 dMat.scale.dc<-double.center.data.matrix(dMat.scale)
 dMat.scale.dc.eigen<-eigen(dMat.scale.dc)
 rownames(dMat.scale.dc.eigen$vectors)<-rownames(dMat)
 #--select eigen-vectors 1 and 2, convert to polar coordinates, and order terms by phi (counter-clockwise)
 dMat.scale.dc.eigen.12<-dMat.scale.dc.eigen$vectors[,1:2]
 dMat.scale.dc.eigen.12.p<-c2p(dMat.scale.dc.eigen.12)
 res<-dMat.scale.dc.eigen.12[order(dMat.scale.dc.eigen.12.p[,"phid"],decreasing=F),]
 return(res)
}

order.terms.using.hclust<-function(goGO)
{
 #--take an IC-subsetted coverage matrix and order terms as follows
 dMat<-make.term.wise.distance.matrix(goGO$covmat)
 dMatHc<-hclust(as.dist(dMat),"single")
 res<-dMatHc$labels[dMatHc$order]
 return(res)
}

tabulate.go.GO.output.aug<-function(goGO)
{
 #--take an output from 'go.GO' and make a table...
 res=NULL
 colNames=NULL
 curres=NULL
 curcolNames=NULL
 goGO.numAnnoInTstSumm<-sapply(goGO$enrichmentRes,FUN=function(x){x$numAnnoInTst})
 goGO.propAnnoInTstSumm<-sapply(goGO$enrichmentRes,FUN=function(x){x$propAnnoInTst})
 goGO.expectedAnnoPropSumm<-sapply(goGO$enrichmentRes,FUN=function(x){x$expectedAnnoProp})
 goGO.adjpSumm<-sapply(goGO$enrichmentRes,FUN=function(x){-log10(p.adjust(x$pval))})

for(curcol in (1:ncol(goGO.numAnnoInTstSumm)))
 {
  curres<-cbind(numAnnoInTest=goGO.numAnnoInTstSumm[,curcol],
                    obsPropAnnoInTst=goGO.propAnnoInTstSumm[,curcol],
                    expPropAnnoInTst=goGO.expectedAnnoPropSumm[,curcol],
                    adjP=goGO.adjpSumm[,curcol])
  curcolNames<-paste(colnames(curres),goGO$catRun[curcol],sep=".")
  res<-cbind(res,curres)
  colNames<-c(colNames,curcolNames)
 }

 colnames(res)<-colNames

 return(res)
}

