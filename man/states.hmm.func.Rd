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
  
  \item{sample}{sample identifier}
  \item{chrom}{chromosome identifier}
  \item{dat}{dataframe with clones in the rows and samples in the
    columns}
  \item{datainfo}{dataframe containing the clones information that is
    used to map each clone of the array to a position on the genome. Has
    to contain columns with names Clone/Chrom/kb containing clone names,
    chromosomal assignment and kb positions respectively}
  \item{vr}{Initial experimental variance}
  \item{maxiter}{Maximum number of iterations}
  \item{aic}{TRUE or FALSE variable indicating whether or nor AIC
    criterion should be used for model selection (see DETAILS)}
  \item{bic}{TRUE or FALSE variable indicating whether or nor BIC
    criterion should be used for model selection (see DETAILS)}
  \item{delta}{numeric vector of penalty factors to use with BIC
    criterion. If BIC is true, delta=1 is always calculated (see DETAILS)}
  \item{nlists}{defaults to 1 when aic=TRUE, otherwise > 1}
  \item{eps}{parameter controlling the convergence of the EM algorithm.}
  \item{print.info}{ print.info = T allows diagnostic information to be
    printed on the screen. }
  \item{diag.prob}{ parameter controlling the construction of the
    initial transition probability matrix. }
}
\seealso{
  \code{\link{aCGH}}
}
\keyword{models}% at least one, from doc/KEYWORDS
