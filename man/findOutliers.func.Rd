\name{findOutliers.func}
\alias{findOutliers.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to identify outlier clones}
\description{
  The function identified the clones that are outliers.
}
\usage{
findOutliers.func(thres = madGenome, factor = 4, statesres = states.bic)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{thres}{Estimate of experimental variability, generally,
    madGenome}
  \item{factor}{Factor indicating how many standard}
  \item{statesres}{The states output of the \code{\link{hmm.run.func}}}
}
\details{
  The outliers are the clones that are dissimilar enough from the clones
  assigned to the same state. Magnitude of the factor determines how
  many MADs away from a median a value needs to be to be declared an
  outlier. Outliers consitent over many samples may indicate
  technological artificat with that clone or possibly copy number
  polymorpism.
}
\value{
  \item{outlier}{Binary matrix with a row for each clone and column for
    each sample. "1" indicates outlier, 0 otherwise.}
  \item{pred.obs.out}{Matrix with a row for each clone and column for
    each sample. The entries are the median value for the state with
    outliers exceluded for all clones but outliers. The value if the
    observed value for the outliers.}
  \item{pred.out}{Matrix with a row for each clone and column for each
    sample. The entries are the median value for the state }
}
\seealso{
  \code{\link{aCGH}}
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH  data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand }
\keyword{models}% at least one, from doc/KEYWORDS
