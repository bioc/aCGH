\name{plotSummaryProfile}
\alias{plotSummaryProfile}
\title{plotSummaryProfile}
\description{
  This function display the genomic events and tests for the differences
  between groups if they are specified.
}
\usage{
plotSummaryProfile(aCGH.obj,
                   response = as.factor(rep("All", ncol(aCGH.obj))),
                   titles = unique(response[!is.na(response)]),
                   X = TRUE, Y = FALSE, maxChrom = 23,
                   chrominfo = human.chrom.info.Jul03)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{an object of aCGH class.}
  \item{response}{
    categorical phenotype to compare. defaults to all the samples being
    analyzed together.
  }
  \item{titles}{
    titles for the groups, defaults to the name of the \code{response}
  }
  \item{X}{
    logical indicating whether X needs to be shown
  }
  \item{Y}{
    logical indicating whether Y needs to be shown
  }  
  \item{maxChrom}{this parameter controls how many chromosomes will be
    plotted, from 1 to maxChrom.}
  \item{chrominfo}{
    a chromosomal information associated with the mapping of the data.
  }
}
\details{
  This function utilizes output of the \code{\link{find.genomic.events}}
  by plotting it and testing between groups. The test are performed
  using kruskal-wallis rank test.
}
% \value{
%   plot containing a number of events
% }
\seealso{
  \code{\link{aCGH}}
  \code{\link{find.genomic.events}}
}
\examples{

data(colorectal)

## Plotting summary of the sample profiles
plotSummaryProfile(colorectal)

}
\keyword{hplot}% at least one, from doc/KEYWORDS
