\name{plotFreqStat}
\alias{plotFreqStat}
\alias{threshold.func}
\alias{changeProp.func}
\alias{changeProp.overall.func}
\alias{table.bac.func}
\alias{lengthGain.na}
\alias{propGain.na}
\alias{lengthLoss.na}
\alias{propLoss.na}
\alias{prop.na}
\alias{gainloss.func}
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
             X = TRUE, Y = FALSE, threshold = TRUE, minChanged = 0,
             rsp.uniq = unique(pheno),
             all = length(rsp.uniq) == 1 && is.null(resT),
             nlim = 1, titles = if (all) "All Samples" else rsp.uniq,
             cutplot = 0, thres = .25, ylm = c(-1, 1), ngrid = 2,
             p.thres = c(.01, .05, .1), mincolors = .1,
             quant.col = .11, numaut = 22, onepage = TRUE,
             colored = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{threshold}{
    \code{threshold} specifies whether significance analysis should be
    performed using original values (F) or thresholded (trimerized)
    values (T).
  }
  \item{minChanged}{
    Only clones that change in at least minChanged proportion of samples
    are plotted; is 0 by default.
  }
  \item{rsp.uniq}{\code{rsp.uniq} specified the codes for the groups of
    interest. Default is the unique levels of the phenotype. Not used
    when \code{all} is T.
  }
  \item{all}{
    \code{all} specifies whether samples should be analyzed by subgroups
    (T) or together (F).
  }
  \item{nlim}{\code{nlim} maximum value for colors.}
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
    threshold. Clone is considered to be gained if it is above the
    threshold and lost if it below negative threshold. Used for plotting
    the gain/loss frequency data as well as for clone screening and for
    significance analysis when \code{threshold} is T.
  }
  \item{ylm}{\code{ylm} vertical limits for the plot}
  \item{ngrid}{
    \code{ngrid} in how many grids red and green colors need to be
    divided.
  }
  \item{p.thres}{
    \code{p.thres} vector of p-value ciut-off to be plotted. computed
    conservatively as the threshold corresponding to a given adjusted
    p-value.
  }
  \item{mincolors}{\code{mincolors} minimum value for colors.}
  \item{quant.col}{
    \code{quant.col} what quantile to use for color display.
  }
  \item{numaut}{\code{numaut} number of the autosomes}
  \item{onepage}{
    \code{onepage} whether all plots are to be plotted on one page or
    different pages.
  }
  \item{colored}{Is plotting in color or not? Default is T.}
}
% \details{
%   ~~ If necessary, more details than the __description__  above ~~
% }
\seealso{
  \code{\link{aCGH}}
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
age.na <- !is.na(age)
colorectal.na <- colorectal[ ,age.na, keep = TRUE ]
dat <- log2.ratios.imputed(colorectal.na)
stat.age <- sapply(1:nrow(dat),
                   function(i) {
                      if (i \%\% 100 == 0)
                         cat(i, "\n")
                      lm.fit <-
                          summary(lm(dat[i,] ~ age[age.na]))
                      c(lm.fit$fstatistic[1],
                        1 - pf(lm.fit$fstatistic[1],
                               lm.fit$fstatistic[2],
                               lm.fit$fstatistic[3])
                       )
                   }
                   )
## Make resT
resT.age <- data.frame(index = 1:ncol(stat.age),
                       teststat = stat.age[ 1, ],
                       rawp = stat.age[ 2, ],
                       adjp = p.adjust(stat.age[ 2, ], "fdr"))
#write out the table of results:
tbl.age <- data.frame(clones.info(colorectal), resT.age)
write.table(tbl.age, "table.age.txt", sep="\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

}
\keyword{htest}
\keyword{hplot}
