\name{colorectal}
\alias{colorectal}
\alias{clones.info.ex}
\alias{log2.ratios.ex}
\alias{pheno.type.ex}
\docType{data}
\title{Colorectal array CGH dataset}
\description{
  The colorectal dataset is an object of class \code{\link{aCGH}}. It
  represents a collection of 124 array CGH profiles of primary
  colorectal tumors and their derived attributes. Each sample was
  measured on the BAC clone DNA microarray with approximate resolution
  of 1.4 Mb per clone. There were approximately 2400 clones spotted on
  the array and each clone was printed in triplicates located
  immediately next to each other. Each array consisted of the 16 (4 by
  4) subarrays. The clones were mapped on the July 03 UCSC freeze. There
  were a number of the discrete and continious phenotypes associated
  with the samples such as age, mutation status for various markers,
  stage, location and so on. All images were quantified and
  normalized by Dr. Taku Tokuyasu using custom image software SPOT and
  postprocessing custom software SPROC.
}
\usage{data(colorectal)}
% \format{
%   The format is the object of class aCGH.
% }
\source{
  These data were generated at Dr. Fred Waldman's lab at UCSF Cancer
  Center by K. Nakao and K. Mehta. The manuscript describing the data
  and the analysis are described in
  High-resolution analysis of DNA copy number alterations in colorectal
  cancer by array-based comparative genomic hybridization,
  \emph{Carcinogenesis}, 2004, Nakao et. al.
}
\references{
  Nakao et. al., High-resolution analysis of DNA copy number
  alterations in colorectal cancer by array-based comparative genomic
  hybridization, \emph{Carcinogenesis}, 2004
  Jain et. al, Fully automatic quantification of microarray image data,
  \emph{Genome Research}, 2003
}
\examples{

data(colorectal)
## WARNING: plotting the heatmap takes some time
plot(colorectal)
plotGenome(colorectal[,1:2])

}
\seealso{
  \code{\link{aCGH}}
  \code{\link{plotGenome}}
}
\keyword{datasets}
