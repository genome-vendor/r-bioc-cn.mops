\name{getReadCountsFromBAM}
\alias{getReadCountsFromBAM}
\title{Calculation of read counts from BAM files.}
\usage{
  getReadCountsFromBAM(BAMFiles, sampleNames, refSeqName,
    WL, mode, parallel = 0)
}
\arguments{
  \item{BAMFiles}{BAMFiles}

  \item{sampleNames}{The corresponding sample names to the
  BAM Files.}

  \item{refSeqName}{Name of the reference sequence that
  should be analyzed. The name must appear in the header of
  the BAM file. If it is not given the function will select
  the first reference sequence that appears in the header
  of the BAM files.}

  \item{WL}{Windowlength. Length of the initial
  segmentation of the genome in basepairs. Should be chosen
  such that on the average 100 reads are contained in each
  segment. If not given, cn.mops will try to find an
  appropiate window length.}

  \item{mode}{Possible values are "paired" and "unpaired",
  whether the mapping algorithm was using a "paired" or
  "unpaired" strategy.}

  \item{parallel}{The number of parallel processes to be
  used for this function. Default=0.}
}
\value{
  An instance of "GRanges", that contains the breakpoints
  of the initial segments and the raw read counts that were
  extracted from the BAM files. This object can be used as
  input for cn.mops and other CNV detection methods.
}
\description{
  Generates the read counts from BAM Files. These counts
  are necessary for CNV detection methods based on depth of
  coverage information. Note that the function is much
  faster, if the BAM files have an index file. The index
  file is assumed to be in the same folder and have an
  identical file name except that ".bai" is appended.

  This function can also be run in a parallel version.
}
\examples{
BAMFiles <- list.files(system.file("extdata", package="cn.mops"),pattern=".bam$",
	full.names=TRUE)
bamDataRanges <- getReadCountsFromBAM(BAMFiles,
					sampleNames=paste("Sample",1:3),WL=5000,mode="unpaired")
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

