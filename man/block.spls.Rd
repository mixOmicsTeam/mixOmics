\name{block.spls}
\encoding{UTF-8}
\alias{block.spls}

\title{N-integration and feature selection with sparse Projection to Latent Structures models (sPLS)}

\description{Integration of multiple data sets measured on the same samples or observations, with variable selection in each data set, ie. N-integration. The method is partly based on Generalised Canonical Correlation Analysis.}


\usage{
block.spls(X,
Y,
indY,
ncomp = 2,
keepX,
keepY,
design,
scheme,
mode,
scale = TRUE,
init ,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
all.outputs = TRUE)
}

\arguments{
\item{X}{A list of data sets (called 'blocks') measured on the same samples. Data in the list should be arranged in matrices, samples x variables, with samples order matching in all data sets.}
\item{Y}{Matrix response for a multivariate regression framework. Data should be continuous variables (see block.splsda for supervised classification and factor reponse)}
\item{indY}{To supply if Y is missing, indicates the position of the matrix response in the list \code{X}}
\item{ncomp}{the number of components to include in the model. Default to 2. Applies to all blocks.}
\item{keepX}{A list of same length as X.  Each entry is the number of variables to select in each of the blocks of X for each component. By default all variables are kept in the model.}
\item{keepY}{Only if Y is provided. Each entry is the number of variables to select in each of the blocks of Y for each component.}
\item{design}{numeric matrix of size (number of blocks in X) x (number of blocks in X) with values between 0 and 1. Each value indicates the strenght of the relationship to be modelled between two blocks; a value of 0 indicates no relationship, 1 is the maximum value. If \code{Y} is provided instead of \code{indY}, the \code{design} matrix is changed to include relationships to \code{Y}. }
\item{scheme}{Either "horst", "factorial" or "centroid". Default = \code{horst}, see reference.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details. Default = \code{regression}.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single"). Default = \code{svd.single}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}
\item{all.outputs}{boolean. Computation can be faster when some specific (and non-essential) outputs are not calculated. Default = \code{TRUE}.}
}


\details{
\code{block.spls} function fits a horizontal sPLS model with a specified number of components per block).
An outcome needs to be provided, either by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
Multi (continuous)response are supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by being disregarded during the cross product computations in the algorithm \code{block.pls} without having to delete rows with missing data. Alternatively, missing data can be imputed prior using the  \code{nipals} function.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References and \code{?pls} for more details).

Note that our method is partly based on sparse Generalised Canonical Correlation Analysis and differs from the MB-PLS approaches proposed by Kowalski et al., 1989, J Chemom 3(1), Westerhuis et al., 1998, J Chemom, 12(5) and sparse variants Li et al., 2012, Bioinformatics 28(19); Karaman et al (2014), Metabolomics, 11(2);
Kawaguchi et al., 2017, Biostatistics.


Variable selection is performed on each component for each block of \code{X}, and for \code{Y} if specified, via input parameter \code{keepX} and \code{keepY}.

Note that if \code{Y} is missing and \code{indY} is provided, then variable selection on \code{Y} is performed by specifying the input parameter directly in \code{keepX} (no \code{keepY} is needed).

}

\value{
\code{block.spls} returns an object of class \code{"block.spls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{indY}{the position of the outcome Y in the output list X.}
\item{ncomp}{the number of components included in the model for each block.}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of each block}
\item{keepY}{Number of variables used to build each component of Y}
\item{variates}{list containing the variates of each block of X.}
\item{loadings}{list containing the estimated loadings for the variates.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each block}

}


\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares.
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.

Tenenhaus A. and Tenenhaus M., (2011), Regularized
Generalized Canonical Correlation Analysis,
Psychometrika, Vol. 76, Nr 2, pp 257-284.

Tenenhaus A., Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V. Variable selection for generalized canonical correlation analysis. \emph{Biostatistics}. kxu001
}

\author{Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao}

\seealso{\code{\link{plotIndiv}}, \code{\link{plotArrow}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.pls}}, \code{\link{block.splsda}} and http://www.mixOmics.org for more details.}


\examples{

# Example with multi omics TCGA study
# -----------------------------
data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna)
# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = c(2)
# set number of variables to select, per component and per data set (this is set arbitrarily)
list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2))
list.keepY = c(rep(10, 2))

TCGA.block.spls = block.spls(X = data, Y = breast.TCGA$data.train$protein,
ncomp = ncomp, keepX = list.keepX, keepY = list.keepY, design = design)
TCGA.block.spls
# in plotindiv we color the samples per breast subtype group but the method is unsupervised!
plotIndiv(TCGA.block.spls, group =  breast.TCGA$data.train$subtype, ind.names = FALSE)
# illustrates coefficient weights in each block
plotLoadings(TCGA.block.spls, ncomp = 1)
plotVar(TCGA.block.spls, style = 'graphics', legend = TRUE)
network(TCGA.block.spls)

}

\keyword{regression}
\keyword{multivariate}
