\name{aCGH.test}
\alias{aCGH.test}
\alias{mt.maxT}
\alias{mt.minP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Testing association of aCGH clones with censored or continuous
  outcomes}
\description{
  \code{aCGH.test} function tests for association of each clone in an
  univariate manner with censored or continous outcome by fitting Cox
  proportional hazards model or linear regression model. There is also
  an alternative to Cox prop. hazards - testing for differences in
  survival curves defined by the groups in the outcome variable using
  the \eqn{G^\rho} family of tests.
}
\usage{
aCGH.test(aCGH.obj, rsp, test = c("survdiff","coxph",
          "linear.regression"), p.adjust.method = "fdr",
           subset = NULL, strt = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{ aCGH object containing clones' log2 ratios. }
  \item{rsp}{ Response variable which is either \code{\link{Surv}}
    object from \code{survival} package or continous outcome.}
  \item{test}{ Currently only three values are allowed - "coxph",
    "survdiff", and "linear.regression", which test for association
    using Cox proportional hazards model, \eqn{G^\rho} family of tests
    (\code{\link{survdiff}}) or linear model.
  }
  \item{p.adjust.method}{ This is a parameter controlling how the
    p-values from the univariate tests are going to be adjusted for
    multiple testing. Default value is Benjamini & Hochberg (1995) FDR
    method. Please refer to \code{\link{p.adjust}} function for more
    help.}
  \item{subset}{Specifies subset index of clones to be tested.}
  \item{strt}{Aptional strata variable for splitting the data in
    different strata.}
  \item{\dots}{ Optional parameters passed further along to each of the
    univariate testing functions.}
}
% \details{
%   ~~ If necessary, more details than the __description__  above ~~
% }
\value{
  A data frame similar to the result returned from
  \code{\link{mt.maxT}} function from \code{multtest} package with
  components:
  \item{index}{
    Vector of row indices, between 1 and \code{nrow(X)}, where rows are
    sorted first according to their adjusted \eqn{p}-values, next their
    unadjusted \eqn{p}-values, and finally their test statistics.
  }
  \item{teststat}{
    Vector of test statistics, ordered according to \code{index}. To get
    the test statistics in the original data order, use
    \code{teststat[order(index)]}.
  }
  \item{rawp}{
    Vector of raw (unadjusted) \eqn{p}-values, ordered according to
    \code{index}.
  }
  \item{adjp}{
    Vector of adjusted \eqn{p}-values, ordered according to
    \code{index}.
  }
}
% \references{ ~put references to the literature/web site here ~ }
\author{Peter Dimitrov}
% \note{ ~~further notes~~ }

%  ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{aCGH}}, \code{\link{Surv}}, \code{\link{mt.maxT}},
  \code{\link{coxph}}, \code{\link{survdiff}}, \code{\link{p.adjust}}}
\examples{
}
\keyword{htest}
