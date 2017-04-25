#--code from Ramnik Xavier's lab at MGH for calculating tissue enrichment scores
#--http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/EnrichmentProfilerDownload/EnrichmentScore.R

require(limma)

#this is the main fucntion
#it takes 2 argument:
#Data, which is the matrix of expression values
#groups,  which is a character vector with group names. Replicates must have the same group name
#and the vector must be in the same order as the columns in Data
GenerateEnrichmentScore<-function(data, groups){
	results <- matrix(NA, nrow=dim(data)[1], ncol=length(unique(groups)))
	rownames(results)<-rownames(data)
	colnames(results)<-unique(groups)
	for (groupX in unique(groups)){
		colX <- RunOneGroupForScore(data, groups, groupX)
		results[,groupX]<-colX
		cat("Group", groupX, "done\n")
		}
	results
	}

#this function will generate the enrichment score for a single
#group.
#Data is the matrix of expression values
#groups is a character vector with group names. Replicates must have the same name
#groupX is the name of the group for which the enrichment score will be calculated
RunOneGroupForScore<-function(data,groups, groupX){
	f<-factor(groups, levels=unique(groups))
	design<-model.matrix(~0+f)
	colnames(design)<-unique(groups)
	mycontrasts<-makeContrasts(contrasts= makeGroupContrastsForOne(groupX,groups), levels=design)
	fit<-lmFit(data, design)
	fit2<-eBayes(contrasts.fit(fit, mycontrasts))
	fit.test<-p.adjust(fit2$p.value, method="bonferroni")
	select <- fit.test > 0.05
	fit.coef <- fit2$coef
	fit.coef[select] <- 0
	apply(fit.coef, 1, sum)
}

#internal function to generate contrasts for limma
makeGroupContrastsForOne<-function(groupX, groups){
	unique.groups<-unique(groups)
	results <- rep("", length(unique.groups)-1)
	counter <- 1
	for (groupY in unique.groups)
		if (groupX!=groupY){
			results[counter]<-paste(groupX,"-", groupY, sep="")
			counter <- counter+1
			}
	results
	}
