\name{impute.HMM}
\alias{impute.HMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Imputing log2 ratios using HMM}
\description{
  Imputing log2 ratios using the output of the HMM segmenttation
}
\usage{
impute.HMM(aCGH.obj, chrominfo = human.chrom.info.Jul03, maxChrom =
23, use.BIC = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{Object of class aCGH.}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{maxChrom}{Highest chromosome to impute.}
  \item{use.BIC}{logical parameter; if true impute missing values based
    on the Hidden Markov Model selected using Bayesian Information
    Criterion impute missing data, otherwise use AIC.}
}
\details{
  See details in \code{\link{aCGH}} discussion.
}
\value{
  Computes and returns the imputed log2 ratio matrix of the aCGH object
  using the output of the Hidden Markov Model segmentation done by
  invoking \code{\link{find.hmm.states}} function.
}
\seealso{
  \code{\link{aCGH}}, \code{\link{find.hmm.states}}, \code{\link{impute.lowess}}.
}
\examples{

datadir <- system.file(package = "aCGH")
datadir <- paste(datadir, "/examples", sep="")

clones.info <-
      read.table(file = file.path(datadir, "clones.info.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
log2.ratios <-
      read.table(file = file.path(datadir, "log2.ratios.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
ex.acgh <- create.aCGH(log2.ratios, clones.info)

## Imputing the log2 ratios 

hmm(ex.acgh) <- find.hmm.states(ex.acgh, aic = TRUE, delta = 1.5)
log2.ratios.imputed(ex.acgh) <- impute.HMM(ex.acgh)

}
\keyword{models}% at least one, from doc/KEYWORDS
