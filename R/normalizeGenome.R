# Copyright (C) 2011 Klambauer Guenter 
# <klambauer@bioinf.jku.at>


#' @title Normalization of NGS data
#' 
#' @description Normalize quantitative NGS data in order to make counts comparable over
#' samples. Scales each samples' reads such that the coverage is even for
#' all samples after normalization. 
#' @param X Matrix of positive real values, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.  Alternatively this can be
#' a GRanges object containing the read counts as values.
#' @param normType normType Type of the normalization technique.
#' Each samples' read counts
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
#' X.norm <- normalizeGenome(X)
#' @return A data matrix of normalized read counts with the same dimensions 
#' as the input matrix X.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export



normalizeGenome <- function(X,normType="poisson",qu=0.25,ploidy){
	YY <- normalizeChromosomes(X,normType=normType,qu=qu,ploidy=ploidy)
	return(YY)
}
