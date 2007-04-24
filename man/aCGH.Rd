\name{aCGH}
\alias{aCGH}
\alias{create.aCGH}
\alias{log2.ratios}
\alias{clones.info}
\alias{is.aCGH}
\alias{dim.aCGH}
\alias{num.clones}
\alias{nrow.aCGH}
\alias{num.samples}
\alias{num.chromosomes}
\alias{ncol.aCGH}
\alias{clone.names}
\alias{row.names.aCGH}
\alias{rownames.aCGH}
\alias{clone.names<-}
\alias{row.names<-.aCGH}
\alias{rownames<-.aCGH}
\alias{col.names.aCGH}
\alias{col.names<-.aCGH}
\alias{colnames.aCGH}
\alias{colnames<-.aCGH}
\alias{sample.names}
\alias{sample.names<-}
\alias{log2.ratios.imputed}
\alias{log2.ratios.imputed<-}
\alias{hmm}
\alias{hmm<-}
\alias{hmm.merged}
\alias{hmm.merged<-}
% \alias{impute.HMM}
\alias{sd.samples}
\alias{sd.samples<-}
\alias{genomic.events}
\alias{genomic.events<-}
\alias{phenotype}
\alias{phenotype<-}
\alias{[.aCGH}
\alias{print.aCGH}
\alias{summary.aCGH}
\alias{plot.aCGH}
\alias{minna}
\alias{maxna}
\alias{corna}
\alias{floor.func}
\alias{lengthNumFunc}
\alias{propNumFunc}
\alias{subset.hmm}
\alias{subset.hmm.merged}
\alias{ex.acgh.hmm}
% \alias{as.eSet}

\title{Class aCGH}
\description{Objects of this class represent batch of arrays of
  Comparative Genomic Hybridization data. In addition to that, there are
  slots for representing phenotype and various genomic events associated
  with aCGH experiments, such as transitions, amplifications,
  aberrations, and whole chromosomal gains and losses.
  Currently objects of class aCGH are represented as S3 classes which
  are named list of lists with functions for accessing elements of that
  list. In the future, it's anticipated that aCGH objects will be
  implemented using S4 classes and methods.
}
%- maybe also 'usage' for other objects documented here.
\details{
  One way of creating objects of class aCGH is to provide the two
  mandatory arguments to \code{create.aCGH} function: \code{log2.ratios}
  and \code{clones.info}. Alternatively aCGH object can be created using
  \code{\link{aCGH.read.Sprocs}} that reads Sproc data files and
  creates object of type aCGH.
}
\note{
  \code{clones.info} slot has to contain a list with at least
  4 columns: Clone (clone name), Target (unique ID, e.g. Well ID), Chrom
  (chromosome number, X chromosome = 23 in human and 20 in mouse, Y
  chromosome = 24 in human and 21 in mouse) and kb (kb position on the
  chromosome).
}
\value{
  
  \item{log2.ratios}{
    Data frame containing the log2 ratios of copy number changes; rows
    correspond to the clones and the columns to the samples (Mandatory).
  }
  \item{clones.info}{
    Data frame containing information about the clones used for
    comparative genomic hybridization. The number of rows of
    \code{clones.info} has to match the number of rows in
    \code{log2.ratios} (Mandatory).
  }
  \item{phenotype}{
    Data frame containing phenotypic information about samples used in
    the experiment generating the data. The number of rows of
    \code{phenotype} has to match the number of columns in
    \code{log2.ratios} (Optional).
  }
  \item{log2.ratios.imputed}{
    Data frame containing the imputed log2 ratios. Calculate this using
    \code{impute.lowess} function; look at the examples below (Optional).
  }
  \item{hmm}{
    The structure of the hmm element is described in
    \code{\link{hmm}}. Calculate this using
    \code{\link{find.hmm.states}} function; look at the examples below
    (Optional).
  }
  \item{hmm}{
    Similar to the structure of the hmm element. Calculate this using
    \code{\link{mergeHmmStates}} function; look at the examples below
    (Optional).
  }
  \item{sd.samples}{
    The structure of the sd.samples element is described in
    \code{\link{computeSD.Samples}}. Calculate this using
    \code{\link{computeSD.Samples}} function; look at the examples
    below (Optional). It is prerequisite that the hmm states are
    estimated first.
  }
  \item{genomic.events}{
    The structure of the genomic.events element is described in
    \code{\link{find.genomic.events}}. Calculate this using
    \code{\link{find.genomic.events}} function; look also at the
    examples below. It is prerequisite that the hmm states and
    sd.samples are computed first. The \code{genomic.events} is used
    widely in variety of plotting functions such as
    \code{\link{plotHmmStates}}, \code{\link{plotFreqStat}}, and
    \code{\link{plotSummaryProfile}}.
  }
  \item{dim.aCGH}{
    returns the dimensions of the aCGH object: number of clones by
    number of samples.
  }
  \item{num.clones}{
    number of clones/number of rows of the log2.ratios data.frame.
  }
  \item{nrow.aCGH}{
    same as \code{num.clones}.
  }
  \item{is.aCGH}{
    tests if its argument is an object of class \code{aCGH}.
  }
  \item{num.samples}{
    number of samples/number of columns of the log2.ratios data.frame.
  }
  \item{nrow.aCGH}{
    same as \code{num.samples}.
  }
  \item{num.chromosomes}{
    number of chromosomes processed and stored in the aCGH object.
  }
  \item{clone.names}{
    returns the names of the clones stored in the clones.info slot of
    the aCGH object.
  }
  \item{row.names.aCGH}{
    same as \code{clone.names}.
  }
  \item{sample.names}{
    returns the names of the samples used to create the aCGH
    object. If the object is created using
    \code{\link{aCGH.read.Sprocs}}, these are the file names of the
    individual arrays.
  }
  \item{col.names.aCGH}{
    same as \code{sample.names}.
  }
  \item{[.aCGH}{
    subsetting function. Works the same way as \code{\link{[.data.frame}}.
  }
%   \item{as.eSet}{
%     function for converting aCGH object to object of \code{\link[Biobase]{eSet-class}}.
%   }
  Most of the functions/slots listed above have assignment operators
  '<-' associated with them.
  
}
\seealso{
  \code{\link{aCGH.read.Sprocs}},
  \code{\link{find.hmm.states}},
  \code{\link{computeSD.Samples}},
  \code{\link{find.genomic.events}},
  \code{\link{plotGenome}},
  \code{\link{plotHmmStates}},
  \code{\link{plotFreqStat}},
  \code{\link{plotSummaryProfile}}
}
\examples{

## Creating aCGH object from log2.ratios and clone info files
## For alternative way look at aCGH.read.Sprocs help

datadir <- system.file(package = "aCGH")
datadir <- paste(datadir, "/examples", sep="")

clones.info <-
      read.table(file = file.path(datadir, "clones.info.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
log2.ratios <-
      read.table(file = file.path(datadir, "log2.ratios.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
pheno.type <-
      read.table(file = file.path(datadir, "pheno.type.ex.txt"),
                 header = TRUE, sep = "\t", quote="", comment.char="")
ex.acgh <- create.aCGH(log2.ratios, clones.info, pheno.type)

## Printing, summary and basic plotting for objects of class aCGH

data(colorectal)
colorectal
summary(colorectal)
sample.names(colorectal)
phenotype(colorectal)
plot(colorectal)

## Subsetting aCGH object

colorectal[1:1000, 1:30]

## Imputing the log2 ratios 

log2.ratios.imputed(ex.acgh) <- impute.lowess(ex.acgh)

## Determining hmm states of the clones
## WARNING: Calculating the states takes some time

##in the interests of time, hmm-finding function is commented out
##instead the states previosuly save are assigned
##hmm(ex.acgh) <- find.hmm.states(ex.acgh)

hmm(ex.acgh) <- ex.acgh.hmm
hmm.merged(ex.acgh) <-
   mergeHmmStates(ex.acgh, model.use = 1, minDiff = .25)

## Calculating the standard deviations for each array

sd.samples(ex.acgh) <- computeSD.Samples(ex.acgh)

## Finding the genomic events associated with each sample

genomic.events(ex.acgh) <- find.genomic.events(ex.acgh)

## Plotting and printing the hmm states

plotHmmStates(ex.acgh, 1)
pdf("hmm.states.temp.pdf")
plotHmmStates(ex.acgh, 1)
dev.off()

## Plotting summary of the sample profiles

plotSummaryProfile(colorectal)

% ## Converting to eSet-class

% as.eSet(colorectal)

}
\author{
  Peter Dimitrov
}
\keyword{classes}% at least one, from doc/KEYWORDS
