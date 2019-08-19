\name{pca}
\encoding{UTF-8}
\alias{pca}

\title{Principal Components Analysis}

\description{
Performs a principal components analysis on the given data matrix that can
contain missing values. If data are complete 'pca' uses Singular Value
Decomposition, if there are some missing values, it uses the NIPALS algorithm.
}

\usage{
pca(X,
ncomp = 2,
center = TRUE,
scale = FALSE,
max.iter = 500,
tol = 1e-09,
logratio = 'none',  # one of ('none','CLR','ILR')
ilr.offset = 0.001,
V = NULL,
multilevel = NULL)
}

\arguments{
\item{X}{a numeric matrix (or data frame) which provides the
data for the principal components analysis. It can contain missing values.}
\item{ncomp}{integer, if data is complete \code{ncomp} decides the number of components and associated
eigenvalues to display from the \code{pcasvd} algorithm and if the data has missing values,
\code{ncomp} gives the number of components to keep to perform the reconstitution of the data using
the NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X), ncol(X))}}
\item{center}{a logical value indicating whether the variables should be shifted to be zero centered.
Alternately, a vector of length equal the number of columns of \code{X} can be supplied.
The value is passed to \code{\link{scale}}.}
\item{scale}{a logical value indicating whether the variables should be scaled to have
unit variance before the analysis takes place. The default is \code{FALSE} for consistency with \code{prcomp}
function, but in general scaling is advisable. Alternatively, a vector of length equal the number of
columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.}
\item{max.iter}{integer, the maximum number of iterations in the NIPALS algorithm.}
\item{tol}{a positive real, the tolerance used in the NIPALS algorithm.}
\item{logratio}{one of ('none','CLR','ILR'). Specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing data. Default to 'none'}
\item{ilr.offset}{When logratio is set to 'ILR', an offset must be input to avoid infinite value after the logratio transform, default to 0.001.}
\item{V}{Matrix used in the logratio transformation id provided.}
\item{multilevel}{sample information for multilevel decomposition for repeated measurements.}
}

\details{
The calculation is done either by a singular value decomposition of the (possibly centered and scaled)
data matrix, if the data is complete or by using the NIPALS algorithm if there is data missing. Unlike
\code{\link{princomp}}, the print method for these objects prints the results in a nice format and the
\code{plot} method produces a bar plot of the percentage of variance explaned by the principal
components (PCs).

When using NIPALS (missing values), we make the assumption that the first (\code{min(ncol(X),} \code{nrow(X)})
principal components will account for 100 \% of the explained variance.

Note that \code{scale= TRUE} cannot be used if there are zero or constant (for \code{center = TRUE}) variables.

According to Filzmoser et al., a ILR log ratio transformation is more appropriate for PCA with compositional data. Both CLR and ILR are valid.

Logratio transform and multilevel analysis are performed sequentially as internal pre-processing step, through \code{\link{logratio.transfo}} and \code{\link{withinVariation}} respectively.

Logratio can only be applied if the data do not contain any 0 value (for count data, we thus advise the normalise raw data with a 1 offset). For ILR transformation and additional offset might be needed.
}

\value{
\code{pca} returns a list with class \code{"pca"} and \code{"prcomp"} containing the following components:
\item{ncomp}{the number of principal components used.}
\item{sdev}{the eigenvalues of the covariance/correlation matrix, though the calculation is actually
done with the singular values of the data matrix or by using NIPALS.}
\item{rotation}{the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors).}
\item{loadings}{same as 'rotation' to keep the mixOmics spirit}
\item{x}{the value of the rotated data (the centred (and scaled if requested) data
multiplied by the rotation/loadings matrix), also called the principal components.}
\item{variates}{same as 'x' to keep the mixOmics spirit}
\item{center, scale}{the centering and scaling used, or \code{FALSE}.}
\item{explained_variance}{explained variance from the multivariate model, used for plotIndiv}
}

\references{
On log ratio transformations:
Filzmoser, P., Hron, K., Reimann, C.: Principal component analysis for compositional data with outliers. Environmetrics 20(6), 621-632 (2009)
Lê Cao K.-A., Costello ME,  Lakis VA, Bartolo, F,Chua XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights into Microbial Communities. PLoS ONE, 11(8): e0160169 (2016).
On multilevel decomposition:
Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde, A.K.: Multivariate paired data analysis: multilevel plsda versus oplsda. Metabolomics 6(1), 119-128 (2010)
Liquet, B., Lê Cao, K.-A., Hocini, H., Thiebaut, R.: A novel approach for biomarker selection and the integration of repeated measures experiments from two assays. BMC bioinformatics 13(1), 325 (2012)
}


\author{Florian Rohart, Kim-Anh Lê Cao, Ignacio González}

\seealso{\code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}},
\code{\link{plotIndiv}}, \code{\link{plotVar}}
and http://www.mixOmics.org for more details.}

\examples{

# example with missing values where NIPALS is applied
# --------------------------------
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

\dontrun{
# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)))

\dontrun{
plotIndiv(pca.res, cex = 0.2,
col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")
}

# variable representation
plotVar(pca.res)

\dontrun{
plotVar(pca.res, rad.in = 0.5, cex = 0.5,style="3d")
}

# example with multilevel decomposition and CLR log ratio transformation (ILR longer to run)
# ----------------

data("diverse.16S")
pca.res = pca(X = diverse.16S$data.TSS, ncomp = 5,
logratio = 'CLR', multilevel = diverse.16S$sample)
plot(pca.res)
plotIndiv(pca.res, ind.names = FALSE, group = diverse.16S$bodysite, title = '16S diverse data',
legend = TRUE)
}
}

\keyword{algebra}
