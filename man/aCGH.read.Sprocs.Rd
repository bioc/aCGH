\name{aCGH.read.Sprocs}
\alias{aCGH.read.Sprocs}
\alias{dotify.names}
\alias{read.Sproc.files}
\alias{extract.clones.info}
\title{Create object of class "aCGH" from Sproc files}
\description{
  This function reads in two-channel Array Complimentary Genomic
  Habrydization Sproc files, flags them for bad quality and missing
  data, and creates object of class \code{aCGH}.
}
\usage{
aCGH.read.Sprocs(fnames, latest.mapping.file = NULL, maxsd = 0.2,
                 minreplic = 2, chrom.remove.threshold = 24,
                 prop.missing = 0.25, sample.names = fnames,
                 sample.quality.threshold = 0.4,
                 cols = c("Log2Rat", "Log2StdDev", "NReplic", "Bad.P")
                 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnames}{a vector of character strings containing the file names
    of each Sproc data file.}
  \item{latest.mapping.file}{The name of an optional file that contains
    newer clone mapping different from the clone mapping used at the
    time when the arrays were created.}
  \item{maxsd}{maximum of standard deviation of log2 ratios used in
    pre-filtering.}
  \item{minreplic}{minimum number of replicates per clone for a single
    chip used to calculate the log2 ratios.}
  \item{chrom.remove.threshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 repsectivelly, in for Mus musculus 20
    and 21, etc.}
  \item{prop.missing}{Clones are screened out and if the proportion
    missing in the samples is \code{prop.missing} they are removed.}
  \item{sample.names}{Sample names. If they are missing, the file names
    are used after stripping the characters after the last dot in the
    filename if one exists; for example 'myfile.txt' becomes myfile.}
  \item{sample.quality.threshold}{Mark those samples that have their
    proportion of missing values \code{sample.quality.threshold}.}
  \item{cols}{
    character vector of length 4 containing the following Sproc file
    column names: log2 ratios, std. deviations of the log2
    ratios, number of replicates for each clone and flags for bad
    clones. Defaults to c("Log2Rat", "Log2StdDev", "NReplic",
    "Bad.P"). Note that all the whitespace characters in the column
    names will be replaced with dots.
  }
  
}
\value{
  Object of class \code{aCGH}.
}
\author{Jane Fridlyand, Peter Dimitrov}
\seealso{
  \code{aCGH}
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

## Testing if creating the object went right. Should all be true.

all(log2.ratios(ex.acgh)[ 1, ] == c(-0.077698 , 0.007389))
clone.name <- "RP11-88I12"
all(log2.ratios(ex.acgh)[ clone.name, ] == c(0.096387 , 0.050999))

}
\keyword{file}% at least one, from doc/KEYWORDS
