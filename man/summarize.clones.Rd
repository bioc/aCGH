\name{summarize.clones}
\alias{summarize.clones}
\title{Extracting summaries for all clones}
\description{
summarize.clones() function is the text equivalent of
\code{\link{plotFreqStat}} function - it summarizes the frequencies of
changes for each clone accross tumors.
}
\usage{
summarize.clones(aCGH.obj, resT = NULL, pheno = rep(1, ncol(aCGH.obj)), rsp.uniq = unique(pheno), minChanged = 0, thres = 0.2, all = length(rsp.uniq) == 1 && is.null(resT), p.thres = c(0.01, 0.05, 0.1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{ ~~Describe \code{aCGH.obj} here~~ }
  \item{resT}{ ~~Describe \code{resT} here~~ }
  \item{pheno}{ ~~Describe \code{pheno} here~~ }
  \item{rsp.uniq}{ ~~Describe \code{rsp.uniq} here~~ }
  \item{minChanged}{ ~~Describe \code{minChanged} here~~ }
  \item{thres}{ ~~Describe \code{thres} here~~ }
  \item{all}{ ~~Describe \code{all} here~~ }
  \item{p.thres}{ ~~Describe \code{p.thres} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plotFreqStat}}, \code{\link{aCGH}} }
\examples{
data(colorectal)
summarize.clones(colorectal)
}
\keyword{models}
