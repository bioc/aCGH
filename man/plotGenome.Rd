\name{plotGenome}
\alias{plotGenome}
\alias{plotGene}
\alias{plotGeneSign}
\alias{maPalette}
\alias{plotCGH.func}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the genome}
\description{
  Basic plot of the log2 ratios for each array ordered along the genome.
}
\usage{
plotGenome(aCGH.obj, samples = 1:num.samples(aCGH.obj), naut = 22,
           Y = TRUE, X = TRUE, data = log2.ratios(aCGH.obj),
           chrominfo = human.chrom.info.Jul03,
           yScale = c(-2, 2), samplenames = sample.names(aCGH.obj),
           ylb = "Log2Ratio")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{an object of class aCGH}
  \item{samples}{vector containing indeces of the samples to be plotted.}
  \item{naut}{number of autosomes in the organism}
  \item{Y}{TRUE if chromosome Y is to be plotted, FALSE otherwise}
  \item{X}{TRUE if chromosome X is to be plotted, FALSE otherwise}
  \item{data}{a matrix containing values to use for plotting. defaults to the \code{log2.ratios(aCGH.obj)}.}
  \item{chrominfo}{a chromosomal information associated with the mapping
    of the data.}
  \item{yScale}{
    Minimum y-scale to use for plotting. Scale is expanded if any of the
    values exceed the positive or negative limit.
  }
  \item{samplenames}{sample names.}
  \item{ylb}{label for the Y-axis.}
}
\seealso{
  \code{\link{aCGH}}
}
\examples{

#plot samples in the order of descending quality 
data(colorectal)
order.quality <- order(sd.samples(colorectal)$madGenome)
postscript("plotGenome.orderByQuality.ps", paper="letter")
par(mfrow=c(2,1))
for(i in order.quality)
   plotGenome(colorectal, samples = i, Y = FALSE)
dev.off()

}
\keyword{hplot}% at least one, from doc/KEYWORDS
