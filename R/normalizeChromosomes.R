# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>

.statmod <- function(x,na.rm=FALSE) {
	if (na.rm){
		z <- table(as.vector(x[!is.na(x)]))
		r <- names(z)[z == max(z)]
		return(as.numeric(r)[1])
	} else {
		if (any(is.na(x))){return(NA)
		} else {
			z <- table(as.vector(x))
			r <- names(z)[z == max(z)]
			return(as.numeric(r)[1])
		}
	}
} 


#' @title Normalization of NGS data. 
#' 
#' Normalize quantitative NGS data in order to make counts comparable over
#' samples, i.e., correcting for different library sizes or coverages. 
#' Scales each samples' reads such that the coverage is even for 
#' all samples after normalization. 
#' 
#' @param X Matrix of positive real values, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region. Alternatively this can be
#' a GRanges object containing the read counts as values.
#' @param chr Character vector that has as many elements as "X" has rows. The
#' vector assigns each genomic segment to a reference sequence (chromosome).
#' @param normType Type of the normalization technique. Each samples'
#' read counts
#' are scaled such that the total number of reads is equal after normalization.
#' By this parameter one can decide to which coverage (i.e. total reads) the 
#' read counts should be normalized. Possible choices are the minimal coverage 
#' ("min"), the mean or median coverage ("mean", "median") or any quantile 
#' ("quant"). If this parameter is set to the value "mode", 
#' the read counts are scaled such that each samples'
#' most frequent value (the "mode") is equal after normalization. If the
#' parameter is set to "poisson" the values are scaled such that the 
#' distribution is (rowwise) close to a Poisson distribution.
#' Possible values are "mean","min","median","quant","poisson, and "mode". 
#' Default = "poisson".
#' @param qu Real value between 0 and 1. Default = 0.25.
#' @param ploidy An integer value for each sample or each column in the read
#' count matrix. At least two samples must have a ploidy of 2. Default = "missing".
#' @examples 
#' data(cn.mops)
#' X.norm <- normalizeChromosomes(X)
#' @return A data matrix of normalized read counts with the same dimensions 
#' as the input matrix X.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


normalizeChromosomes <-
		function(X,chr,normType="poisson",qu=0.25,ploidy){
	if (!(normType %in% c("mean","min","median","quant","mode","poisson"))){
		stop(paste("Set TO of normalization to \"mean\"",
						"\"min\", \"median\", \"quant\" or \"mode\"."))
	}
	input <- X
	returnGRanges <- FALSE
	if(class(X)=="GRanges"){
		returnGRanges <- TRUE
		X <- IRanges::as.matrix(IRanges::values(X)) 	
	}	
	if (is.vector(X)){X <- matrix(X,nrow=1)}
	
	if (missing(chr)){
		chr <- rep("undef",nrow(X))
	}
	if (missing(ploidy)){
		ploidy <- rep(2,ncol(X))
	}
	if (any(ploidy!=as.integer(ploidy))){
		stop("Ploidy values must be integers!")
	}
	if (length(ploidy)!=ncol(X)){
		stop("Length of the ploidy vector does not match the number of", 
				"columns of the read count matrix!")
	}
	ploidy <- as.integer(ploidy)
	if (!length(which(ploidy>=2))){
		stop("At least two diploid samples must be contained in the data.")
	}
	ploidy[ploidy==0] <- 0.05
	Xorig <- X
	
	# Sequencing data matrix
	# vector of chromosome - length equal to rows of X
	if (length(chr)!=nrow(X)){
		stop("Length of \"chr\" must be equal to number of rows of \"X\".")}
	chr <- (as.character(chr))
	
	YY <- matrix(0,nrow=nrow(Xorig),ncol=ncol(Xorig))
	
	ploidy2flag <- FALSE
	ploidy2median <- c()
	
	for (pp in unique(c(2,ploidy))){
		X <- Xorig[,ploidy==pp,drop=FALSE]
		if (ncol(X)==1){
			Y <- X
			
		} else {
			
			Y <- matrix(0,nrow=nrow(X),ncol=ncol(X))
			for (l in (unique(chr))){
				chrIdx <- which(chr==l)
				Ytmp <- X[chrIdx, ,drop=FALSE]
				idxSG <- apply(Ytmp,1,function(x) all(x<1))
				Ytmp[idxSG, ] <- NA
				
				if (nrow(Ytmp) > 1){
					
					mappedreads <- colSums(Ytmp,na.rm=TRUE)
					if (any(mappedreads==0)){
						warning(paste("There exists a reference sequence with zero reads"
										,"for some samples."))
						mappedreads <- pmax(mappedreads, 1)
					}
					
					if (normType == "min"){
						correctwiththis <-  min(mappedreads,na.rm=TRUE)/mappedreads
						
					} else if (normType=="mean"){
						correctwiththis <-  mean(mappedreads,na.rm=TRUE)/mappedreads
					} else if (normType=="median"){
						cat("normalizing to median \n")
						correctwiththis <-  median(mappedreads,na.rm=TRUE)/mappedreads
					} else if (normType=="quant"){
						correctwiththis <-  quantile(mappedreads,probs=qu,
								na.rm=TRUE)/mappedreads
					} else if (normType=="mode"){
						ssm <- apply(Ytmp,2,function(x) .statmod(x[x!=0],na.rm=TRUE) )
						asm <- .statmod(Ytmp[Ytmp!=0],na.rm=TRUE)
						correctwiththis <- asm/ssm
						
					} else if (normType=="poisson"){
						#browser()
						correctwiththis <-  median(mappedreads,na.rm=TRUE)/mappedreads
						YYtmp <- t(t(Ytmp)*correctwiththis)
						v2m <- apply(YYtmp,1,var)/rowMeans(YYtmp)
						uut <- quantile(v2m,probs=0.95,na.rm=TRUE)
						dd <- density(v2m,na.rm=TRUE,from=0,
								to=uut,n=uut*100)
						#mv2m <- median(v2m,na.rm=TRUE)
						mv2m <- dd$x[which.max(dd$y)]
						if (is.finite(mv2m)) {correctwiththis <- correctwiththis*1/mv2m}
					} else {
						stop("normType not known.")
					}
					if (any(!is.finite(correctwiththis))){
						warning(paste("Normalization for reference sequence ",l,"not", 
										"applicable, because at least one sample has zero",
										"reads."))
						correctwiththis <- rep(1,ncol(X))
					}
					Ytmp <- t(t(Ytmp)*correctwiththis)
					Ytmp[idxSG, ] <- 0
				}
				
				Y[chrIdx, ] <- Ytmp
				
				
			} # over chr
			
			if (!ploidy2flag){
				ploidy2flag <- TRUE
				ploidy2median <- median(Y[!idxSG, ],na.rm=TRUE)
			}
			
		}
		#browser()
		
		YY[,ploidy==pp] <- Y*ploidy2median/median(Y[!idxSG, ],na.rm=TRUE)*pp/2
	}
	
	rownames(YY) <- rownames(Xorig)
	colnames(YY) <- colnames(Xorig)
	
	if (returnGRanges){
		values(input) <- YY
		return(input)
	} else {
		return(YY)
	}
}

