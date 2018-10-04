\name{spca}
\encoding{UTF-8}
\alias{spca}

\title{Sparse Principal Components Analysis}

\description{
Performs a sparse principal components analysis to perform variable selection by using singular value decomposition.}

\usage{
spca(X, ncomp = 2, center = TRUE, scale = TRUE,
keepX = rep(ncol(X),ncomp), max.iter = 500,
tol = 1e-06, logratio = 'none',  # one of ('none','CLR')
multilevel = NULL)
}

\arguments{
\item{X}{a numeric matrix (or data frame) which provides the
data for the sparse principal components analysis.}
\item{ncomp}{integer, the number of components to keep.}
\item{center}{a logical value indicating whether the variables should be shifted to be zero centered.
Alternatively, a vector of length equal the number of columns of \code{X} can be supplied.
The value is passed to \code{\link{scale}}.}
\item{scale}{a logical value indicating whether the variables should be scaled to have
unit variance before the analysis takes place. The default is \code{TRUE}. See details.}
\item{max.iter}{integer, the maximum number of iterations to check convergence in each component.}
\item{tol}{a positive real, the tolerance used in the iterative algorithm.}
\item{keepX}{numeric vector of length ncomp, the number of variables to keep in loading vectors. By default all
variables are kept in the model. See details.}
\item{logratio}{one of ('none','CLR'). Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing data. Default to 'none'}
\item{multilevel}{sample information for multilevel decomposition for repeated measurements.}

}

\details{
The calculation employs singular value decomposition of the (centered and scaled) data matrix and LASSO to generate
sparsity on the loading vectors.

\code{scale= TRUE} is highly recommended as it will help obtaining orthogonal sparse loading vectors.

\code{keepX} is the number of variables to keep in loading vectors. The difference between number of columns
of \code{X} and \code{keepX} is the degree of sparsity, which refers to the number of zeros in each loading vector.

Note that \code{spca} does not apply to the data matrix with missing values. The biplot function for \code{spca} is not available.

According to Filzmoser et al., a ILR log ratio transformation is more appropriate for PCA with compositional data. Both CLR and ILR are valid.

Logratio transform and multilevel analysis are performed sequentially as internal pre-processing step, through \code{\link{logratio.transfo}} and \code{\link{withinVariation}} respectively.

Logratio can only be applied if the data do not contain any 0 value (for count data, we thus advise the normalise raw data with a 1 offset). For ILR transformation and additional offset might be needed.
}

\value{
\code{spca} returns a list with class \code{"spca"} containing the following components:
\item{ncomp}{the number of components to keep in the calculation.}
\item{varX}{the adjusted cumulative percentage of variances explained.}
\item{keepX}{the number of variables kept in each loading vector.}
\item{iter}{the number of iterations needed to reach convergence for each component.}
\item{rotation}{the matrix containing the sparse loading vectors.}
\item{x}{the matrix containing the principal components.}
}

\references{
Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via regularized
low rank matrix approximation. \emph{Journal of Multivariate Analysis} \bold{99}, 1015-1034.
}

\author{Kim-Anh Lê Cao, Fangzhou Yao, Leigh Coonan}

\seealso{\code{\link{pca}} and http://www.mixOmics.org for more details.}

\examples{
data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
spca.rat

## variable representation
plotVar(spca.rat, cex = 0.5)
plotVar(spca.rat,style="3d")

\dontrun{
## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3],
group = as.numeric(liver.toxicity$treatment[, 3]))
plotIndiv(spca.rat, cex = 0.01,
col = as.numeric(liver.toxicity$treatment[, 3]),style="3d")


# example with multilevel decomposition and CLR log ratio transformation
# ----------------

data("diverse.16S")
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 5,
logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, group = diverse.16S$bodysite, title = '16S diverse data',
legend=TRUE)

}
}

\keyword{algebra}
