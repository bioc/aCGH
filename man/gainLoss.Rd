\name{gainLoss}
\alias{gainLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute proportion of gains and losses for each clones}
\description{
 This function outputs lists containing proportion of gains and losses for each clone. 
}
\usage{
gainLoss(dat, cols, thres=0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{log2ratios of the relevant array CGH object}
  \item{cols}{indeces of the samples to use}
  \item{thres}{global or tumor-specific threshold. defaults to 0.25}
  }
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  \item{gainP}{Vector of proportion gained for each clones}
  \item{lossP}{Vector of proportion lost for each clones}
 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jane Fridlyand}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plotFreqStat}}}
\examples{

data(colorectal)

## Use mt.maxT function from multtest package to test
## differences in group means for each clone grouped by sex
##use only clones with show gain or loss in at least 10\% of the samples
colnames(phenotype(colorectal))
sex <- phenotype(colorectal)$sex
sex.na <- !is.na(sex)
colorectal.na <- colorectal[ ,sex.na, keep = TRUE ]
factor <- 2.5
minChanged <- 0.1
gainloss <- gainLoss(log2.ratios(colorectal.na), cols=1:ncol(colorectal.na), thres=factor*sd.samples(colorectal.na)$madGenome)
ind.clones.use <- which(gainloss$gainP >= minChanged | gainloss$lossP>= minChanged)
#create filtered dataset
colorectal.na <- colorectal.na[ind.clones.use,keep=TRUE]
dat <- log2.ratios.imputed(colorectal.na)
resT.sex <- mt.maxT(dat, sex[sex.na],test = "t.equalvar", B = 1000)


## Plot the result along the genome
plotFreqStat(colorectal.na, resT.sex, sex[sex.na],factor=factor,titles = c("Male", "Female"))


}
\keyword{htest}

