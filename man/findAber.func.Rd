\name{findAber.func}
\alias{findAber.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to  determines focal aberrations}
\description{
  The function identifies clones that are focal aberrations.
}
\usage{
findAber.func(maxClones = 1, maxLen = 1000, statesres = states.bic)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{maxClones}{Maximum number of clones assigned to the same state
    which can be considered to be focal aberrations}
  \item{maxLen}{Maximum lengeth of the region containing clones assigned
    to the state so that those clones can be considered to be focal
    aberrations}
  \item{statesres}{The states output of the \code{\link{hmm.run.func}}}
}
\details{
  The focal aberrations are the one or more clones assigned to the
  state different from the states of the surrounding clones.  They may
  indicate copy number polymorphisms or interesting high or low focal
  changes.
}
\value{
  \item{aber}{Binary matrix with a row for each clone and column for
    each sample. 1 indicates presence of a focal aberrations, 0 lack
    of such.}
}
\references{"Application of Hidden Markov Models to the analysis of the array CGH  data", Fridlyand et.al., JMVA, 2004
}
\author{Jane Fridlyand}
\keyword{models}% at least one, from doc/KEYWORDS
