\name{findTrans.func}
\alias{findTrans.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Funtion identifying the transitions}
\description{
  This function identifies the start and end of the states (regions with
  the constant estimated copy number).
}
\usage{
findTrans.func(outliers = res1$outliers, aber = res2$aber, statesres =
states.bic)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{outliers}{Binary matrix of the outliers (generally output of the
    \code{\link{findOutliers.func}})}
  \item{aber}{Binary matrix of the focal aberrations (generally output
    of the \code{\link{findAber.func}})}
  \item{statesres}{The states output of the \code{\link{hmm.run.func}}}
}
\details{
  The transitions end is placed at the last non-focal aberration clone
  of the contiguous region containing clones belonging to the same state
  and transitions start is placed at the first non-focal aberration
  clone of the contiguous region containing clones belonging to the same
  state.
}
\value{
  \item{trans.matrix}{Matrix with a row for each clone and column for
    each sample. The starts of the states are indicated by "1" , the end
    are by "2" and the focal aberrations are coded as "3"}
  \item{translen.matrix}{Matrix with a row for each clone and column for
    each sample. The entries are the length of the region to which a
    clone belongs. Zero length is assigned to the focal
    aberrations. This output may be buggy at the moment.}
}
\seealso{
  \code{\link{aCGH}}
}
\references{
  Application of Hidden Markov Models to the analysis of the array CGH
  data, Fridlyand et.al., \emph{JMVA}, 2004.
}
\author{Jane Fridlyand}
\keyword{models}% at least one, from doc/KEYWORDS
