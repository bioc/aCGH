\name{plotHmmStates}
\alias{plotHmmStates}
\alias{plotHmmStatesNew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plotting the estimated hmm states and log2 ratios for each sample.
}
\description{
  This function displays the estimated hmm states and log2 ratios for
  each sample.
}
\usage{
plotHmmStates(aCGH.obj, sample.ind, chr = 1:num.chromosomes(aCGH.obj),statesres =
              hmm(aCGH.obj)$states.hmm[[1]],
              maxChrom = 23, chrominfo = human.chrom.info.Jul03,
              yScale = c(-2,2), samplenames = sample.names(aCGH.obj))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{object of class aCGH}
  \item{sample.ind}{index of the sample to be plotted relative to the
    data matrix (i.e. column index in the file)}
  \item{statesres}{matrix containing states informations. defaults to the states selected using the first  model selection criterionof \code{aCGH.obj}}
  \item{chr}{vector of chromosomes to be plotted}
  \item{yScale}{specified scale for Y-exis}
  \item{maxChrom}{highest chromosome to show}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{samplenames}{vector of sample names}
}
\details{
  Each chromosome is plotted on a separate  page and contains two
  figures. The top figure shows the observed log2ratios and the bottom
  figure shows predicted values for all clones but outliers which show
  observed values. The genomic events are indicated on both figures as
  following. The first clone after transition is indicated with solid
  blue line and the last clone after transitions is shown with dotted
  green line. Focal aberrations clones are colored orange,
  amplifications are colored red and outliers are yellow.
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand}
\seealso{
  \code{\link{aCGH}}
  \code{\link{find.hmm.states}}
  \code{\link{plotGenome}}
}
\examples{

data(colorectal)
plotHmmStates(colorectal, 1)

}
\keyword{hplot}% at least one, from doc/KEYWORDS
