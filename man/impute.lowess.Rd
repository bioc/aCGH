\name{impute.lowess}
\alias{impute.lowess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Imputing log2 ratios}
\description{
  Imputing log2 ratios
}
\usage{
impute.lowess(aCGH.obj, chrominfo = human.chrom.info.Jul03, maxChrom =
23, smooth = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{Object of class aCGH.}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{maxChrom}{Highest chromosome to impute.}
  \item{smooth}{smoothing parameter for the lowess procedure}
}
\details{
  There are two main reasons to impute data. One is that given that imputation is reasonable, one can increase the analytical power and improve results. Another, more practical, is that at the moment many widely used fuctions in R do not support missing values. While procedures such as kNN imputations is widely used for gene expression data, it is more powerful to take advantage of the genomic structure of the array CGH data and use a smoother. Note that we perform only one pass os smoothing. If there still remain missing values, they are imputed by the median on the chromosome or chromosomal arm where applicable,
}
\value{
  Computes and returns the imputed log2 ratio matrix of the aCGH object.
}
\seealso{
  \code{\link{aCGH}}
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

log2.ratios.imputed(ex.acgh) <- impute.lowess(ex.acgh)

}
\keyword{models}% at least one, from doc/KEYWORDS
