\name{plotvalGenome.func}
\alias{plotvalGenome.func}
\alias{plotvalChrom.func}
\alias{plotValGenome}
\alias{plotValChrom}
\alias{plotChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{clustering and heatmap}
\description{
  This function clusters samples and shows their heatmap
}
\usage{
plotvalGenome.func(aCGH.obj,
                   response = as.factor(rep("All", ncol(aCGH.obj))),
                   chrominfo = human.chrom.info.Jul03, cutoff=1,
                   lowCol = "red", highCol = "green", midCol = "black",
                   ncolors = 50, byclass = FALSE, showaber = FALSE,
                   amplif = 1, homdel = -0.75,
                   samplenames = sample.names(aCGH.obj),
                   vecchrom = 1:23, titles = "Image Plot",
                   methodS = "ward", dendPlot = TRUE, imp = TRUE,
                   categoricalPheno = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aCGH.obj}{object of class \code{aCGH} here}
  \item{response}{phenotype of interest. defaults to the same phenotype
    assigned to all samples}
  \item{chrominfo}{a chromosomal information associated with the mapping
    of the data}
  \item{cutoff}{maximum absolute value. all the values are floored to
    +/-cutoff depending on whether they are positive of
    negative. defaults to 1}
  \item{ncolors}{number of colors in the grid. input to
    \code{\link{maPalette}}. defaults to 50}
  \item{lowCol}{color for the low (negative) values. input to
    \code{\link{maPalette}}. defaults to "red"}
  \item{highCol}{color for the high (positive) values. input to
    \code{\link{maPalette}}. defaults to "green"}
  \item{midCol}{color for the values close to 0. input to
    \code{\link{maPalette}}. defaults to "black"}
  \item{byclass}{logical indicating whether samples should be clustered within each level of the phenotype or overall. defaults to F}
  \item{showaber}{logical indicating whether high level amplifications and homozygous deletions should be indicated on the plot. defaults to F}
  \item{amplif}{positive value that all observations equal or exceeding it  are marked by yellow dots indicating high-level changes. defaults to 1}
  \item{homdel}{negative value that all observations equal or below it  are marked by light blue dots indicating homozygous deletions. defaults to -0.75}
  \item{samplenames}{sample names}
  \item{vecchrom}{vector of chromosomal indeces to use for clustering and to display. defaults to 1:23}
  \item{titles}{plot title. defaults to "Image Plots"}
  \item{methodS}{clustering method to cluster samples. defaults to "ward"}
  \item{dendPlot}{logical indicating whether dendogram needs to be
    drawn. defaults to T.}
  \item{imp}{logical indicating whether imputed or original values should be used. defaults to T, i.e. imputed.}
  \item{categoricalPheno}{logical indicating whether phenotype is
    categorical. Continious phenotypes are treated as "no groups" except
    that their values are dispalyed.defaults to TRUE.
  }
}
\details{
  This functions is a more flexible version of the
  \code{\link{heatmap}}. It can cluster within levels of categorical
  phenotype as well as all of the samples while displaying phenotype
  levels in different colors. It also uses any combination of
  chromosomes that is requested and clusters samples based on these
  chromosomes only. It draws the chromosomal boundaries and displays
  high level changes and homozygous deletions. If phenotype if not
  categical, its values may still be displayed but groups are not formed
  and byclass = F.
  Image plot has the samples reordered according to clustering order.
}
% \value{
% }
%\references{ ~put references to the literature/web site here ~ }
\seealso{
  \code{\link{aCGH}}
  \code{\link{heatmap}}
}
\examples{ 
data(colorectal)

#cluster all samples using imputed data on all chromosomes (autosomes and X):

plotvalGenome.func(colorectal)

#cluster samples within sex groups based on 3 chromosomes individually. 
#use non-imputed data and  do not show dendogram. Indicate amplifications and 
#homozygous deletions.

plotvalGenome.func(colorectal, response = phenotype(colorectal)$sex,
                   byclass = TRUE, showaber = TRUE, vecchrom = c(4,8,9),
                   dendPlot = FALSE, imp = FALSE)

#cluster samples based on each chromosome individualy and display age. Show
#gains in red and losses in green. Show aberrations and use values < -1
#to identify homozgous deletions. Do not show dendogram.

pdf("plotimages.pdf", width = 11, height = 8.5)
for (i in 1:23)
    plotvalGenome.func(colorectal,
                       response = phenotype(colorectal)$age,
                       chrominfo = human.chrom.info.Jul03,
                       cutoff = 1, ncolors = 50, lowCol="green",
                       highCol="red", midCol="black", byclass = FALSE,
                       showaber = TRUE, homdel = -1, vecchrom = i,
                       titles = "Image Plot", methodS = "ward",
                       dendPlot = FALSE, categoricalPheno = FALSE)
dev.off()

}
\keyword{hplot}
\keyword{cluster}
