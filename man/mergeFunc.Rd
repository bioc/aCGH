\name{mergeFunc}
\alias{mergeFunc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Funtion to merge states based on their state means}
\description{
  This function takes an output of hmm.run.fun (states.hmm) with a
  particular model selection criterion andietratively merges the states
  which means are closer than a supplied threshold.
}
\usage{
mergeFunc(statesres = states.bic, minDiff = 0.1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{statesres}{the sublist of the states.hmm list output from
    find.hmm.states for a given model selection crterion}
  \item{minDiff}{The states whose predicted values are less than minDiff
    apart are merged into one state and all the predicited values are
    recomputed}
}
\details{
  This function is intended to reduce effect of the possible small
  magnitude technological artifcats on the structure determination.
}
\value{
  List containing states.hmm object is returned.
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand}
\keyword{models}% at least one, from doc/KEYWORDS
