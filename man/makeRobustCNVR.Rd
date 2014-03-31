\name{makeRobustCNVR}
\alias{makeRobustCNVR}
\title{Calculates robust CNV regions.}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}

  \item{robust}{Robustness parameter. The higher the value,
  the more samples are required to have a CNV that confirms
  the CNV region. Setting this parameter to 0 restores the
  original CNV regions. (Default=0.5)}

  \item{minWidth}{The minimum length measured in genomic
  regions a CNV region has to span in order to be called. A
  parameter of the segmentation algorithm. (Default=4).}
}
\value{
  \code{makeRobustCNVR} returns a "CNVDetectionResult"
  object containing new values in the slot "cnvr".
}
\description{
  This generic function calculates robust CNV regions by
  segmenting the I/NI call per genomic region.
}
\details{
  cn.mops usually reports a CNV region if at least one
  individual has a CNV in this region. For some
  applications it is useful to find more common CNV
  regions, i.e., regions in which more than one sample has
  a CNV. The I/NI call measures both signal strength and
  how many sample show an abnormal copy number, therefore
  segmentation of the I/NI call can provide robust CNV
  regions.
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
rr <- calcIntegerCopyNumbers(makeRobustCNVR(r,robust=0.1,minWidth=3))
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

