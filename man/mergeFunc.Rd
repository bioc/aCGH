\name{mergeFunc}
\alias{mergeFunc}
\alias{mergeHmmStates}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Funtion to merge states based on their state means}
\description{
  \code{mergeFunc} takes the output of hmm.run.func (or find.hmm.states)
  with a particular model selection criterion and iteratively merges the
  states with means closer than a supplied threshold.
  \code{mergeHmmStates} is a frontend for \code{mergeFunc} using
  \code{\link{aCGH}} object.
}
\usage{
mergeHmmStates(aCGH.obj, model.use = 1, minDiff = 0.25)
mergeFunc(statesres, minDiff = 0.1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{
    Object of class \code{\link{aCGH}}.
  }
  \item{statesres}{the sublist of the states.hmm list output from
    \code{\link{find.hmm.states}} for a given model selection crterion
  }
  \item{minDiff}{The states whose predicted values are less than minDiff
    apart are merged into one state and all the predicited values are
    recomputed.
  }
  \item{model.use}{
    Model selection criterion to use, See \code{\link{find.hmm.states}}.
  }
}
\details{
  This function is intended to reduce effect of the possible small
  magnitude technological artifacts on the structure determination.
}
\value{
  List containing states.hmm object is returned.
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand}
\seealso{\code{\link{aCGH}}, \code{\link{find.hmm.states}}}

\keyword{models}% at least one, from doc/KEYWORDS
