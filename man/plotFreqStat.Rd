\name{plotFreqStat}
\alias{plotFreqStat}
%- \alias{threshold.func}
\alias{changeProp.func}
\alias{changeProp.overall.func}
\alias{table.bac.func}
\alias{lengthGain.na}
\alias{propGain.na}
\alias{lengthLoss.na}
\alias{propLoss.na}
\alias{prop.na}
\alias{create.resT}
\alias{plotFreqStatColors}
\alias{plotFreqStatGrey}
\alias{plotFreqGivenStat}
\alias{plotfreqGivenStatFinalColors}
\alias{plotfreq.stat.final.func}
\alias{plotfreq.stat.chrom.final.func}
\alias{plotfreq.givenstat.final.colors.func}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{frequency plots and significance analysis}
\description{
  The main application of this function is to plot the frequency of
  changes.
}
\usage{
plotFreqStat(aCGH.obj, resT = NULL, pheno = rep(1, ncol(aCGH.obj)),
             chrominfo = human.chrom.info.Jul03,
             X = TRUE, Y = FALSE,
             rsp.uniq = unique(pheno),
             all = length(rsp.uniq) == 1 && is.null(resT),
             titles = if (all) "All Samples" else rsp.uniq,
             cutplot = 0, thres = .25, factor = 2.5, ylm = c(-1, 1),
             p.thres = c(.01, .05, .1), numaut = 22, onepage = TRUE,
             colored = TRUE)
}
\arguments{
  \item{aCGH.obj}{Object of class \code{aCGH}}
  \item{resT}{Data frame having the same structure as the result of
    applying \code{\link{mt.maxT}} or \code{\link{mt.minP}} functions
    from Bioconductor's \code{multtest} package for multiple testing.
    The result is a data frame including the following 4 components:
    'index', 'teststat', 'rawp' and 'adjp'.
  }
  \item{pheno}{phenotype to compare.}
  \item{chrominfo}{
    Chromosomal information. Defaults to
    \code{\link{human.chrom.info.Jul03}}
  }
  \item{X}{Include X chromosome? Defaults to yes.}
  \item{Y}{Include Y chromosome? Defaults to no.}
  \item{rsp.uniq}{\code{rsp.uniq} specified the codes for the groups of
    interest. Default is the unique levels of the phenotype. Not used
    when \code{all} is T.
  }
  \item{all}{
    \code{all} specifies whether samples should be analyzed by subgroups
    (T) or together (F).
  }
   \item{titles}{
    \code{titles} names of the groups to be used. Default is the unique
    levels of the \code{pheno}.
  }
  \item{cutplot}{only clones with at least \code{cutplot} frequency of
    gain and loss are plotted.
  }
  \item{thres}{\code{thres} is either a vector providing unique
    threshold for each sample or a vector of the same length as number
    of samples (columns in \code{data}) providing sample-specific
    threshold. If \code{aCGH.obj} has non-null sd.samples, then \code{thres} is automatically replaced by \code{factor} times madGenome of \code{aCGH} object. Clone is considered to be gained if it is above the
    threshold and lost if it below negative threshold. Used for plotting
    the gain/loss frequency data as well as for clone screening and for
    significance analysis when \code{threshold} is TRUE.Defaults to 0.25
  }
  \item{factor}{\code{factor} specifies the number by which experimental variability should be multiplied. used only when sd.samples(\code{aCGH.obj}) is NULL. Defaults to 2.5}
  \item{ylm}{\code{ylm} vertical limits for the plot}
  \item{p.thres}{
    \code{p.thres} vector of p-value ciut-off to be plotted. computed
    conservatively as the threshold corresponding to a given adjusted
    p-value.
  }
  \item{numaut}{\code{numaut} number of the autosomes}
  \item{onepage}{
    \code{onepage} whether all plots are to be plotted on one page or
    different pages. When more than 2 groups are compared, we recommend multiple pages.
  }
  \item{colored}{Is plotting in color or not? Default is TRUE.}
}
\examples{

data(colorectal)

## Use mt.maxT function from multtest package to test
## differences in group means for each clone grouped by sex
colnames(phenotype(colorectal))
sex <- phenotype(colorectal)$sex
sex.na <- !is.na(sex)
colorectal.na <- colorectal[ ,sex.na, keep = TRUE ]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na], test = "t", B = 1000)

## Plot the result along the genome
plotFreqStat(colorectal.na, resT.sex, sex[sex.na],
             titles = c("Male", "Female"))

## Adjust the p.values from previous exercise with "fdr"
## method and plot them
resT.sex.fdr <- resT.sex
resT.sex.fdr$adjp <- p.adjust(resT.sex.fdr$rawp, "fdr")
plotFreqStat(colorectal.na, resT.sex.fdr, sex[sex.na],
             titles = c("Male", "Female"))

## Derive statistics and p-values for testing the linear association of
## age with the log2 ratios of each clone along the samples

age <- phenotype(colorectal)$age
age.na <- which(!is.na(age))
age <- age[age.na]
colorectal.na <- colorectal[, age.na]
stat.age <- aCGH.test(colorectal.na, age, test = "linear.regression", p.adjust.method = "fdr")

#separate into two groups: < 70 and > 70 and plot freqeuncies of gain and loss
#for each clone. Note that statistic plotted corresponds to linear coefficient
#for age variable

plotFreqStat(colorectal.na, stat.age, ifelse(age < 70, 0, 1), titles =
             c("Young", "Old"), X = FALSE, Y = FALSE)
}
\keyword{htest}
\keyword{hplot}
