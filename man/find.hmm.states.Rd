\name{find.hmm.states}
\alias{find.hmm.states}
\alias{hmm.run.func}
\alias{mu1.func}
\alias{as.matrix.dist}
\alias{plotChrom.hmm.func}
\alias{plotCGH.hmm.func}
\alias{smoothData.func}
\alias{thresholdData.func}
\alias{plotChrom.samples.func}
\alias{plotChrom.grey.samples.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Determines states of the clones}
\description{
  This function runs unsupervised HMM algorithm and produces the
  essentual state information which is used for the subsequent structure
  determination.
}
\usage{
hmm.run.func(dat, datainfo = clones.info, vr = 0.01, maxiter = 100,
             aic = TRUE, bic = TRUE, delta = NA)
find.hmm.states(aCGH.obj, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{
    object of class \code{\link{aCGH}}.
  }
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
  \item{...}{
    All the parameters that can be passed to hmm.run.func except dat and
    datainfo.
  }
}
\details{
  One or more model selection criterion is used to determine number of
  states on each chromosomes. If several are specified, then a separate
  matrix is produced for each criterion used. Delta is a fudge factor in
  BIC criterion: $\delta BIC(\gamma) = \log RSS(\gamma) +
  q_{\gamma}\delta\log n/n.$ Note that delta = NA leads to conventional
  BIC. (Broman KW, Speed TP (2002) A model selection approach for the
  identification of quantitative trait loci in experimental crosses
  (with discussion). J Roy Stat Soc B 64:641-656, 731-775 )

  find.hmm.states(aCGH.obj, ...) uses aCGH object instead of log2 ratios
  matrix dat. Equivalent representation (assuming normally distributed
  residuals) is to write -loglik(gamma) = n/2*log(RSS)(gamma) and then
  bic= -loglik+log(n)*k*delta/2 and aic = -loglik+2*k/2
}
\value{
  Two lists of lists are returned. Each list contains information on the
  states with each of the specified model selection criteria. E.g., if
  AIC = T, BIC = T and delta = c(1.5), then each list will contain
  three lists corresponding to AIC, BIC(1) and BIC(1.5) as the 1st,2nd
  and 3rd lists repsectively. If AIC is used, it always comes first
  followed by BIC and then deltaBIC in the order of delta vector.
  \item{states.hmm}{Each of the sublists contains 2+ 6*n columns where
    the first two columns contain chromosome and kb positions for each
    clone in the dataset supplied followed up by 6 columns for each
    sample where n = number of samples.

    column 1 = state

    column 2 = smoothed value for a clone

    column 3 = probability of being in a state

    column 4 = predicted value of a state

    column 5 = dispersion

    column 6 = observed value
  }
  \item{nstates.hmm}{Each of the sublists contains a matrix with each
    row corresponding to a chromosome and each column to a sample. The
    entries indicate how many different states were identified for a
    given sample on a given chromosome}
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004 }
\author{Jane Fridlyand}
\section{WARNING}{When algortihm fails to fit an HMM for a given number
  of states on a chromosome, it prints a warning.}
\seealso{
  \code{\link{aCGH}}
}
\examples{

datadir <- system.file("examples", package = "aCGH")
latest.mapping.file <-
      file.path(datadir, "human.clones.info.Jul03.csv")
ex.acgh <-
    aCGH.read.Sprocs(dir(path = datadir,
                     pattern = paste("*", "txt", sep = "\."),
                     full.names = TRUE), latest.mapping.file,
                     chrom.remove.threshold = 23)
ex.acgh

hmm(ex.acgh) <- find.hmm.states(ex.acgh, aic = TRUE, delta = 1.5)
summary(ex.acgh)

}
\keyword{models}
