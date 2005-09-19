\name{states.hmm.func}
\alias{states.hmm.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{A function to fit unsupervised Hidden Markov model}
\description{
  This function is a workhorse of \code{\link{find.hmm.states}}. It
  operates on the individual chromosomes/samples and is not called
  directly by users.
}
\usage{
states.hmm.func(sample, chrom, dat, datainfo = clones.info, vr = 0.01,
                maxiter = 100, aic = FALSE, bic = TRUE, delta = 1,
                nlists = 1, eps = .01, print.info = FALSE,
                diag.prob = .99)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{sample}{ ~~Describe \code{sample} here~~ }
  \item{chrom}{ ~~Describe \code{chrom} here~~ }
  \item{dat}{ ~~Describe \code{dat} here~~ }
  \item{datainfo}{ ~~Describe \code{datainfo} here~~ }
  \item{vr}{ ~~Describe \code{vr} here~~ }
  \item{maxiter}{ ~~Describe \code{maxiter} here~~ }
  \item{aic}{ ~~Describe \code{aic} here~~ }
  \item{bic}{ ~~Describe \code{bic} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{nlists}{ ~~Describe \code{nlists} here~~ }
  \item{eps}{ parameter controlling the convergence of the EM algorithm }
  \item{print.info}{ print.info = T allows diagnostic information to be
    printed on the screen. }
  \item{diag.prob}{ parameter controlling the construction of the
    initial transition probability matrix. }
}
\seealso{
  \code{\link{aCGH}}
}
\keyword{models}% at least one, from doc/KEYWORDS
