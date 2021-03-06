\name{getSegmentReadCountsFromBAM}
\alias{getSegmentReadCountsFromBAM}
\title{Calculation of read counts from BAM files for predefined segments.}
\usage{
  getSegmentReadCountsFromBAM(BAMFiles, GR, sampleNames,
    mode, parallel = 0, BAIFiles = NULL)
}
\arguments{
  \item{BAMFiles}{BAMFiles}

  \item{sampleNames}{The corresponding sample names to the
  BAM Files.}

  \item{GR}{A genomic ranges object that contains the
  genomic coordinates of the segments.}

  \item{mode}{Possible values are "paired" and "unpaired",
  whether the mapping algorithm was using a "paired" or
  "unpaired" strategy.}

  \item{parallel}{The number of parallel processes to be
  used for this function. Default=0.}

  \item{BAIFiles}{The names of the BAI files that belong to
  the BAM files. The vector has to be in the same order as
  the vector BAMFiles. If the BAI files have the same name
  as the BAM files, only with ".bai" attached, this
  parameter needs not be set. (Default = NULL).}
}
\value{
  An instance of "GRanges", that contains the breakpoints
  of the initial segments and the raw read counts that were
  extracted from the BAM files. This object can be used as
  input for cn.mops and other CNV detection methods.
}
\description{
  Generates the read counts from BAM Files for predefined
  segments. This is the appropiate choice for exome
  sequencing data, where the bait regions, target regions
  or exons are the predefined segments. These counts are
  necessary for CNV detection methods based on depth of
  coverage information. Note that the function is much
  faster, if the BAM files have an index file. The index
  file is assumed to be in the same folder and have an
  identical file name except that ".bai" is appended.

  This function can also be run in a parallel version.
}
\examples{
BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
	full.names=TRUE)
gr <- GRanges(c("20","20"),IRanges(c(60000,70000),c(70000,80000)))
bamDataRanges <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr,mode="unpaired")
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

