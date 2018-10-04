\name{plotDiablo}
\encoding{UTF-8}
\alias{plotDiablo}

\title{Graphical output for the DIABLO framework}

\description{Function to visualise correlation between components from different data sets}

\usage{
plotDiablo(x, ncomp = 1, legend = TRUE, legend.ncol, ...)
}

\arguments{
\item{x}{object of class inheriting from \code{"block.splsda"}.}
\item{ncomp}{Which component to plot calculated from each data set. Has to be lower than the minimum of \code{object$ncomp}}
\item{legend}{boolean. Whether the legend should be added. Default is TRUE.}
\item{legend.ncol}{Number of columns for the legend. Default to \code{min(5,nlevels(x$Y))}}
\item{\dots}{not used}
}


\details{
The function uses a plot.data.frame to plot the component \code{ncomp} calculated from each data set to visualise whether DIABLO (block.splsda) is successful at maximising the correlation between each data sets' component. The lower triangular panel indicated the Pearson's correlation coefficient, the upper triangular panel the scatter plot.
}




\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery. Submitted.
}

\author{Amrit Singh}

\seealso{\code{\link{block.splsda}} and http://www.mixOmics.org/mixDIABLO for more details.}

\examples{

data('breast.TCGA')
Y = breast.TCGA$data.train$subtype

data = list(mrna =  breast.TCGA$data.train$mrna,
mirna =  breast.TCGA$data.train$mirna, prot =  breast.TCGA$data.train$protein)

# set number of component per data set
ncomp = 3
# set number of variables to select, per component and per data set (arbitrarily set)
list.keepX = list(mrna = rep(20, 3), mirna = rep(10,3), prot = rep(10,3))

# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design

BC.diablo = block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
plotDiablo(BC.diablo, ncomp = 1)
}

\keyword{regression}
\keyword{multivariate}
\value{none}
