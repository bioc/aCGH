\name{summarize.clones}
\alias{summarize.clones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracting summary information for all clones}
\description{
 \code{summarize.clones} function is the text equivalent of
\code{\link{plotFreqStat}} function - it summarizes the frequencies of
changes for each clone accross tumors and when available assigns statistics. The resulting table can be easily exported.
}
\usage{
summarize.clones(aCGH.obj, resT = NULL, pheno = rep(1, ncol(aCGH.obj)), rsp.uniq = unique(pheno), thres = 0.25, factor = 2.5, all = length(rsp.uniq) == 1 && is.null(resT), titles = if (all) "all" else rsp.uniq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{\code{aCGH.obj} object here}
  \item{resT}{Data frame having the same structure as the result of
    applying \code{\link{mt.maxT}} or \code{\link{mt.minP}} functions
    from Bioconductor's \code{multtest} package for multiple testing.
    The result is a data frame including the following 4 components:
    'index', 'teststat', 'rawp' and 'adjp'.Default is the unique levels of the phenotype. Not used when \code{all} is TRUE.}
  \item{pheno}{phenotype to compare}
  \item{rsp.uniq}{\code{rsp.uniq} specified the codes for the groups of
    interest. Default is the unique levels of the phenotype. Not used
    when \code{all} is TRUE.}
  \item{thres}{\code{thres} is either a vector providing unique
    threshold for each sample or a vector of the same length as number
    of samples (columns in \code{data}) providing sample-specific
    threshold. If \code{aCGH.obj} has non-null \code{sd.samples}, then threshold is automatically replaced by tumor-specific sd multiplied by \code{factor}. Clone is considered to be gained if it is above the
    threshold and lost if it below negative threshold. Defaults to 0.25}
  \item{factor}{\code{factor} specifies the number by which experimental variability should be multiples. used only when tumor specific variability in \code{aCGH.obj} is not NULL. Defaults to 2.5}
  \item{all}{\code{all} specifies whether samples should be analyzed by subgroups   (TRUE) or together (FALSE)}
  \item{titles}{\code{titles} names of the groups to be used. Default is the unique levels of the \code{pheno}.}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  Returns matrix containg the following information for each clones: annotation (same as in \code{clones.info}), number and proportion of samples where clone is present,gained and lost; and the same in each group if more than one group. Additionally, if significance comparison has been done, value of the statistic, unadjusted p-value and adjusted p-values are included for each clone.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jane Fridlyand}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plotFreqStat}}, \code{\link{aCGH}}}
\examples{
data(colorectal)
summarize.clones(colorectal)
}
\keyword{models}
