\name{human.chrom.info.Jul03}
\alias{human.chrom.info.Jul03}
\alias{ex.acgh.hmm}
\docType{data}
\title{
  Basic Chromosomal Information for UCSC Human Genome Assembly July 2003
  freeze
}
\description{
  This dataset contains basic chromosomal information for UCSC Human
  Genome Assembly July 2003 freeze.
}
\usage{
human.chrom.info.Jul03 is loaded automatically with
\code{\link{aCGH}} package.
}
\format{
  A data frame with 24 observations on the following 3 variables.
  \describe{
    \item{chrom}{Chromosomal index, X is coded as 23 and Y as 24.}
    \item{length}{Length of each chromosome in kilobases.}
    \item{centromere}{
      Location of the centromere on the chromosome (kb).
    }
  }
}
\details{
  This file is used for many plotting functions and needs to correspond
  to \code{\link{clones.info}} mapping file. The centromeric location is
  approximately extimated by taking mid-point between the last
  fish-mapped clone on the p-arm and the first fish-mapped clone on the
  q-arm using relevant UCSC freeze.
  For an alternative freeze, one needs to manually create a 3-column
  file of the format described above.
}
\source{
  http://genome.ucsc.edu/cgi-bin/hgText
}
% \references{
  
% }
\keyword{datasets}
