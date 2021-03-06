\name{iniCall}
\alias{iniCall}
\title{This generic function returns the informative/non-informative call of
 a CNV detection method stored in an instance of
\code{\link{CNVDetectionResult-class}}.
 The I/NI call is a measure for a genomic
segment across all samples, whether this segment is a CNV region
(informative) or a normal genomic region (non-informative).}
\arguments{
  \item{object}{An instance of "CNVDetectionResult"}
}
\value{
  \code{iniCall} returns a "GRanges" object containing the
  individual calls.
}
\description{
  This generic function returns the
  informative/non-informative call of a CNV detection
  method stored in an instance of
  \code{\link{CNVDetectionResult-class}}.  The I/NI call is
  a measure for a genomic segment across all samples,
  whether this segment is a CNV region (informative) or a
  normal genomic region (non-informative).
}
\examples{
data(cn.mops)
r <- cn.mops(X[1:100,1:5])
iniCall(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at}
}

