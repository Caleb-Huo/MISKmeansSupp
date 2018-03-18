##' generate S for MISKmeans
##'
##' generate S for MISKmeans
##' @title generateS
##' @param seed random seed
##' @param S number of studies
##' @param Types number of omics types. I.e. Gene expression, DNA methylation
##' @param k number of clusters
##' @param meanSamplesPerK mean samples per cluster
##' @param nModule number of modules. A module is a group of genes.
##' @param meanGenesPerModule number of genes per module
##' @param Gmean gene expression template follows N(Gmean,Gsd^2)
##' @param Gsd gene expression template follows N(Gmean,Gsd^2)
##' @param sigma1 noise 1
##' @param sigma2 noise 2
##' @param sigma3 noise 3
##' @param G0 number of noise genes
##' @param nconfounder number of confounders
##' @param nrModule number of modules for confounding variables
##' @param rMeanSubtypes number of subtypes defined by confounding variables
##' @param diffmu effect size difference for subtype predictive genes
##' @param fold how to vary subtype predictive gene signal. 1: original. 0: no signal.
##' @param rho para for inverse Wishart distribution.
##' @param df.prior para for inverse Wishart distribution.
##' @param groupProb subtype predictive genes have prior group information. By prob 1-groupProb, the information will be altered.
##' @return alist
##' @author Caleb
##' @export
##' @examples
##' Sdata =generateS(seed=15213,S=3,Types=2, k=3)
##' #
##' length(Sdata)
##' dim(Sdata[[1]]$d)
##' sum(Sdata[[1]]$subPredictGeneUnion)
generateS <- function(seed=15213,S=3,Types=2,k=3,meanSamplesPerK=c(40,40,30),nModule=30,meanGenesPerModule=30,
                      Gmean=9,Gsd=2,sigma1=1,sigma2=1,sigma3=1,G0=5000,
                      nconfounder=4,nrModule=20,rMeanSubtypes=3,diffmu=1,fold=rep(1,S),
                      rho = 0.5,df.prior = 100,groupProb=1){

  set.seed(seed)

  ## main function:
  ## generate S study
  finalResult <- NULL
  result=NULL

  ## number of correlated genes for the subtype predictive genes
  nGClust=NULL

  ## gene template for the confounding part
  ## number of correlated genes for the confounding part
  rnGClust = NULL

  ## determine number of subtypes for the confounding variable
  rk <- rep(k,nconfounder)

  ## prepare samples for the subtypes
  nall = meanSamplesPerK
  cumnall = cumsum(nall)
  n = sum(nall)
  Sindex = NULL
  label = numeric(n)
  for(i in 1:length(nall)){
    if(i==1){
      Sindex[[i]] = 1:cumnall[i]
    } else {
      Sindex[[i]] = (cumnall[i-1]+1):cumnall[i]
    }
    ## label the samples
    label[Sindex[[i]]] = i
  }

  muLayer1_list <- replicate(Types,list())
  for(s in 1:S){
	  result_s <- NULL
	  for(atype in 1:Types){
		  
		if(length(muLayer1_list[[atype]]) == 0){
			muLayer1 <- NULL
		} else {
			muLayer1 <- muLayer1_list[[atype]]
		}		  
		  
	    dataMain = generateStructure(nModule,Gmean,Gsd,Sindex,sigma1,sigma2,rho = rho,
			df.prior = df.prior,meanGenesPerModule=meanGenesPerModule,
			amuLayer1=muLayer1,anGClust=nGClust,constrain=TRUE,diffmu=diffmu,fold=fold[s])

	    resData = dataMain$data
	    subtypeModuleIndex = dataMain$subtypeModuleIndex
		muLayer1_list[[atype]] <- dataMain$amuLayer1
	    #gplots::heatmap.2(resData, col="greenred" ,trace="none",Rowv=TRUE,  Colv=NA,keysize=1.3)

	    ## next step, simulate correlated random genes, with rk~POI(k) clusters
	    ## the biology plausibility is that there are other factors would influnce the gene
	    ## expression, e.g. gender, race, geological information, disease grade, level of hormone
	    ## then exp~UNIF(), for each of the rk subclass
	    ## generate m unrelated subtypes, for each m, there should be mmoduls

	    if(length(rk)!=0){
		    for (i in 1:length(rk)){
		      rSindex = rPartationSample(n,rk[i])
		  	  rdatai = generateStructure(nrModule,Gmean,Gsd,rSindex,sigma1,sigma2,rho = rho,df.prior = df.prior,meanGenesPerModule=meanGenesPerModule)
		      resData = rbind(resData,rdatai$data)
		    }
	    }
	    dim(resData)
		confounderIndex <- (nrow(dataMain$data) + 1):nrow(resData)

	    #gplots::heatmap.2(resData, col="greenred" ,trace="none",Rowv=TRUE,  Colv=NA,keysize=1.3)

	    ## finally, add some house keeping genes
	    ## sigma3 = 1, error of the random genes
	    templateRandom = rnorm(G0,Gmean,Gsd)

	    randomPart  = matrix(NA,nrow=G0,ncol=n)  ## 20*n
	    for(i in 1:G0){
	      randomPart[i,] = rnorm(n=n,mean=templateRandom[i],sd=sigma3)
	    }

	    resData = rbind(resData,randomPart)

		randomIndex <- (nrow(resData) - nrow(randomPart) + 1):nrow(resData)
		
		result_s[[atype]] <- list(S=resData,label=label,subtypeModuleIndex=subtypeModuleIndex,confounderIndex=confounderIndex,randomIndex=randomIndex)
	    result[[s]] <- result_s
	  }	## end of loops for different omics types  	  	  
  } ## end of loops for S studies

  ## preparing grouping information  
  group <- NULL  
  for(atype in 1:length(result_s)){
	  aresult <- result_s[[atype]]
	  agroup <- aresult$subtypeModuleIndex
	  noiseIndex <- c(aresult$confounderIndex, aresult$randomIndex)
	  for(g in 1:length(agroup)){
		  bgroup <- agroup[[g]]
		  relinkIndex <- rbinom(length(bgroup),1,groupProb) == 0
		  bgroup[relinkIndex] = sample(noiseIndex, sum(relinkIndex), replace=FALSE)
		  agroup[[g]] <- bgroup
	  }
	  group[[atype]] <- agroup
  }


  ## generate omics data input for MISKmeans
  for(s in 1:S){
    result_s <- result[[s]]
    omics <- lapply(result_s,function(x) x$S)
    labels <- result_s[[1]]$label
    subPredictGene <- lapply(result_s,function(x) unlist(x$subtypeModuleIndex))
    subPredictGeneUnion <- NULL
    for(atype in 1:length(result_s)){
  	  ap <- nrow(omics[[atype]])
  	  if(atype==1){
  		subPredictGeneUnion <- 1:ap %in% subPredictGene[[atype]]
  		groupUnion <- group[[atype]]
  		d <- t(omics[[atype]])
  	  } else {
  	  	subPredictGeneUnion <- c(subPredictGeneUnion, 1:ap %in% subPredictGene[[atype]])
  		agroup <- group[[atype]]
  		for(g in 1:length(agroup)){
  			groupUnion[[g]] <- c(groupUnion[[g]], (agroup[[g]] + ap))
  		}
  		d <- cbind(d,t(omics[[atype]]))
  	  }
    } ## end of loops for atype
	finalResult[[s]] <- list(omics=omics, d=d, labels=labels, group=group, subPredictGene=subPredictGene, groupUnion=groupUnion, subPredictGeneUnion=subPredictGeneUnion)
  }## end of loops for S studies
  

  return(finalResult)
}
