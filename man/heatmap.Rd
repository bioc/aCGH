\name{heatmap}
\alias{heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default method for plotting array CGH objects}
\description{
  Clusters samples and produces heatmapp of the observed log2ratios.
}
\usage{
heatmap(x, imp = TRUE, Rowv = NA, Colv = NULL, distfun = dist,
        hclustfun = hclust, add.expr, symm = FALSE,
        revC = identical(Colv, "Rowv"), scale = "none",
        na.rm = TRUE, margins = c(5, 5), ColSideColors,
        RowSideColors, cexRow = 0.2 + 1 / log10(nr),
        cexCol = 0.2 + 1 / log10(nc), labRow = NULL,
        labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
        verbose = getOption("verbose"), methodR = "ward",
        methodC = "ward", zlm = c(-0.5, 0.5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of the \code{\link{aCGH}} object}
  \item{imp}{logical variable indicating whether
    \code{log2.ratios.imputed} or \code{log2.ratios} slot of
    \code{\link{aCGH}} should be used. Defaults to imputed value (TRUE).
  }
  \item{Rowv}{determines if and how the row dendrogram should be computed
    and reordered.  Either a 'dendrogram' or a vector of values
    used to reorder the row dendrogram or 'NA' to suppress any
    row dendrogram (and reordering) or by default, 'NULL'
  }
  \item{Colv}{determines if and how the column dendrogram should be
    reordered.  Has the same options as the \code{Rowv} argument above
    and additionally when \code{x} is a square matrix, \code{Colv = "Rowv"}
    means that columns should be treated identically to the rows.
  }
  \item{distfun}{ function used to compute the distance (dissimilarity) between
    both rows and columns. Defaults to 'dist'.}
  \item{hclustfun}{function used to compute the hierarchical clustering when
    'Rowv' or 'Colv' are not dendrograms.  Defaults to 'hclust'}
  \item{add.expr}{expression that will be evaluated after the call to 'image'.
    Can be used to add components to the plot.}
  \item{symm}{logical indicating if 'x' should be treated *symm*etrically;
    can only be true when 'x' is a square matrix.}
  \item{revC}{logical indicating if the column order should be 'rev'ersed
    for plotting, such that e.g., for the symmetric case, the
    symmetry axis is as usual.}
  \item{scale}{character indicating if the values should be centered and
    scaled in either the row direction or the column direction,
    or none.  The default is "row" if \code{symm} false, and
    "none" otherwise.}
  \item{na.rm}{logical indicating whether 'NA''s should be removed.}
  \item{margins}{numeric vector of length 2 containing the margins (see
    'par(mar= *)') for column and row names, respectively.}
  \item{ColSideColors}{(optional) character vector of length 'ncol(x)'
    containing the color names for a horizontal side bar that may
    be used to annotate the columns of 'x'.}
  \item{RowSideColors}{(optional) character vector of length 'nrow(x)'
    containing the color names for a vertical side bar that may
    be used to annotate the rows of 'x'.}
  \item{cexRow, cexCol}{positive numbers, used as 'cex.axis' in for the row or
    column axis labeling.  The defaults currently only use number
    of rows or columns, respectively.}
  \item{labRow, labCol}{character vectors with row and column labels to use;
    these default to 'rownames(x)' or 'colnames(x)',
    respectively.}
  \item{main, xlab, ylab}{main, x- and y-axis titles;}
  \item{verbose}{logical indicating if information should be printed.}
  \item{methodR}{method to use for clustering rows. defaults to "ward"}
  \item{methodC}{method to use for clustering columns. defaults to "ward"}
  \item{zlm}{all the values greater or equal  than zlm are set to zlm -
    0.01. a;; value less or equal to -zlm are set to -zlm + 0.01
  }
  \item{\dots}{additional arguments passed on to 'image', e.g., 'col'
    specifying the colors.
  }
}
\details{
  This function is almost identical to the \code{\link{heatmap}} in base
  R. The slight modifications are that (1) a user can specify clustering
  method for rows and columns; (2) all the values outside specified
  limits are floored to be 0.01 less than a limit; (3) default values
  are different. Note that using default option of \code{imp} (TRUE)
  produces nicer looking plots as all missing values are removed.
}
\value{
  Invisibly, a list with components
  \item{crowInd}{row index permutation vector as returned by
    \code{\link{order.dendrogram}}
  }
  \item{colInd}{row index permutation vector as returned by
    \code{\link{order.dendrogram}}
  }
}
\references{heatmap function in base R}
\seealso{
  \code{\link{aCGH}}
  \code{\link{plotvalGenome.func}}
}
\examples{

#default plotting method for the aCGH object
data(colorectal)
plot(colorectal)

#to produce smoother looking heatmap, use imp = T: this will use imputed
#slot of aCGH object

plot(colorectal, imp = TRUE)

}
\keyword{hplot}
\keyword{cluster}
