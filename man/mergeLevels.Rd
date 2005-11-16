\name{mergeLevels}
\alias{mergeLevels}
\alias{combine.func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mergeLevels}
\description{
  Merging of predicted levels for array CGH data and similar.
}
\usage{
mergeLevels(vecObs,vecPred,pv.thres=0.0001,ansari.sign=0.05,thresMin=0.05,thresMax=0.5,verbose=1,scale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vecObs}{ Vector of observed values, i.e. observed log2-ratios }
  \item{vecPred}{ Vector of predicted values, i.e. mean or median of levels predicted by segmentation algorithm }
  \item{pv.thres}{ Significance threshold for Wilcoxon test for level merging }
  \item{ansari.sign}{ Significance threshold for Ansari-Bradley test }
\item{thresMin}{ merge if segment medians are closer than thresMin , defaiult is 0.05}
\item{thresMax}{ don't merge if segment medians are further than thresMax (unless needs to be merged for a different reason: wilcoxon test), default is .5}
\item{verbose}{ if 1, progress is printed}
\item{scale}{ whether thresholds are on the log2ratio scale and thus need to be converted to the copy number. default is TRUE
}
}
\details{
  mergeLevels takes a vector of observed log2-ratios and predicted log2ratios and merges levels that are not significantly distinct.
}
\value{
  \item{vecMerged }{Vector with merged values. One merged value returned for each predicted/observed value }
  \item{mnNow }{Merged level medians}
  \item{sq }{Vector of thresholds, the function has searched through to find optimum. Note, these thresholds are based on copy number transformed values}
  \item{ansari }{The p-values for the ansari-bradley tests for each threshold in sq}

}
\references{ Willenbrock H, Fridlyand J. (2005). A comparison study: applying segmentation to array CGH data for downstream analyses.
Bioinformatics. 2005 Sep 14; [Epub ahead of print]  }
\author{ Hanni Willenbrock (\email{Hanni@cbs.dtu.dk}) and Jane Fridlyand (\email{jfridlyand@cc.ucsf.edu}) }
\note{ vecObs and vecPred must have same length and observed and predicted value for a given probe should have same position in vecObs and vedPred. The function assumes that log2-ratios are supplied }

\examples{
# Example data of observed and predicted log2-ratios
vecObs <- c(rep(0,40),rep(0.6,15),rep(0,10),rep(-0.4,20),rep(0,15))+rnorm(100,sd=0.2)
vecPred <- c(rep(median(vecObs[1:40]),40),rep(median(vecObs[41:55]),15),
  rep(median(vecObs[56:65]),10),rep(median(vecObs[66:85]),20),rep(median(vecObs[86:100]),15))

# Plot observed values (black) and predicted values (red)
plot(vecObs,pch=20)
points(vecPred,col="red",pch=20)

# Run merge function
merge.obj <- mergeLevels(vecObs,vecPred)

# Add merged values to plot
points(merge.obj$vecMerged,col="blue",pch=20)

# Examine optimum threshold
merge.obj$sq
}
\keyword{htest}

