\name{aCGH.process}
\alias{aCGH.process}
\title{Process data in aCGH object}
\description{
  This function takes object of class \code{aCGH}, and filters clones based on their mapping information and proportion missing. It also average duplicated clones and reports quality statistic.
}
\usage{
aCGH.process(aCGH.obj, chrom.remove.threshold = 24,
                 prop.missing = 0.25, sample.quality.threshold = 0.4,
                 unmapScreen=TRUE, dupRemove = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{Object of class aCGH}
  \item{chrom.remove.threshold}{Chromosomes are ordered and numbered as
    usual, except for X and Y chromosome, which in for Homo sapiens
    genome have numbers 23 and 24 repsectivelly, in for Mus musculus 20
    and 21, etc.}
  \item{prop.missing}{Clones are screened out and if the proportion
    missing in the samples is \code{prop.missing} they are removed.}
  \item{sample.quality.threshold}{Mark those samples that have their
    proportion of missing values \code{sample.quality.threshold}.}
  \item{unmapScreen}{Indicator for whether clones with incomplete mapping information should be removed from the dataset. Note that leaving them in may cause plotting routines fail. Defaults to TRUE}
   \item{dupRemove}{Indicator for whether clones with duplicate names should be averaged and removed from the dataset leaving only one occurence of each duplicated set.Defaults to TRUE}
	

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


datadir <- system.file(package = "aCGH")
datadir <- paste(datadir, "/examples", sep="")

clones.info <-
      read.table(file = file.path(datadir, "clones.info.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
log2.ratios <-
      read.table(file = file.path(datadir, "log2.ratios.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
pheno.type <- read.table(file = file.path(datadir, "pheno.type.ex.txt"),header = TRUE, sep = "\t", quote="", comment.char="")
ex.acgh <- create.aCGH(log2.ratios, clones.info, pheno.type)
ex.acgh <-
    aCGH.process(ex.acgh, chrom.remove.threshold = 23, prop.missing = .25, 
sample.quality.threshold = .4, unmapScreen=TRUE, dupRemove = FALSE)
ex.acgh

}
\keyword{file}% at least one, from doc/KEYWORDS
