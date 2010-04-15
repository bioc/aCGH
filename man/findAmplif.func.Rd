\name{findAmplif.func}
\alias{findAmplif.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Function to determine high level amplifications}
\description{
  This function identifies high level amplifications by considering the
  height, the width of an amplicon relative to the urrounding
  clones. Only narrow peaks much higher than its neigbors are considered
  as high level amplifications.
}
\usage{
findAmplif.func(absValSingle = 1, absValRegion = 1.5, diffVal1 = 1,
diffVal2 = 0.5, maxSize = 10000, translen.matr, trans.matr, aber,
outliers, pred, pred.obs, statesres)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{absValSingle}{A clone is declared to be an amplification if it
    is a focal aberration or an outlier and its value exceeds
    absValSingle}
  \item{absValRegion}{A clone is an amplification if if a clone belong
    to a region with width less than maxSize and  observed value for a
    clones is greater than absValRegion}
  \item{diffVal1}{Clone is an amplification if it is an aberration and
    greater by diffVal1 than max of the two surrounding stretches}
  \item{diffVal2}{Clone is an amplification if it is an outlier, its
    observed values is greater by diffVal2 than max of the two
    surrounding stretches}
  \item{maxSize}{The clones  may not be declared as amplifications if
    they belong to the states with spanning more than maxSize}
  \item{translen.matr}{State length matrix. The output of the
    \code{\link{findTrans.func}}}
  \item{trans.matr}{Transition matrix. The output of the
    \code{\link{findTrans.func}}}
  \item{aber}{Aberration matrix. The output of the
    \code{\link{findAber.func}}}
  \item{outliers}{Outliers matrix. The output of the
    \code{\link{findOutliers.func}}}
  \item{pred}{Predicted values matrix. The output of the
    \code{\link{findOutliers.func}}}
  \item{pred.obs}{Predicted values matrix with observed values assigned
    to the outliers. The output of the \code{\link{findOutliers.func}}}
  \item{statesres}{The states output of the \code{\link{hmm.run.func}}}
}
\details{
  Note that all the distances are in Megabases and all the heights are
  on log2ratio scale.
}
\value{
  \item{amplif.matrix}{Binary matrix with a row for each clone and
    column for each sample. "1" indicates amplification}
  ...
}
\seealso{
  \code{\link{aCGH}}
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH  data, Fridlyand et.al., \emph{JMVA}, 2004}
\author{Jane Fridlyand}
\keyword{models}% at least one, from doc/KEYWORDS
