\encoding{UTF-8}
\name{plot.tune}
\alias{plot.tune.block.splsda}
\alias{plot.tune.splsda}
\alias{plot.tune}

\title{Plot for model performance}

\description{
Function to plot performance criteria, such as classification error rate or balanced error rate on a tune.splsda result.
}

\usage{
\method{plot}{tune.splsda}(x, optimal = TRUE, sd = TRUE, col, \dots)

\method{plot}{tune.block.splsda}(x, sd = TRUE, col, \dots)
}


\arguments{
\item{x}{an \code{tune.splsda} object.}
\item{optimal}{If TRUE, highlights the optimal keepX per component}
\item{sd}{If 'nrepeat' was used in the call to 'tune.splsda', error bar shows the standard deviation if sd=TRUE}
\item{col}{character (or symbol) color to be used, possibly vector. One color per component.}
\item{\dots}{Further arguments sent to \code{\link{xyplot}} function.}
}

\details{
\code{plot.tune.splsda} plots the classification error rate or the balanced error rate from x$error.rate, for each component of the model.
A lozenge highlights the optimal number of variables on each component.

\code{plot.tune.block.splsda} plots the classification error rate or the balanced error rate from x$error.rate, for each component of the model. The error rate is ordered by increasing value, the yaxis shows the optimal combination of keepX at the top (e.g. `keepX on block 1'_`keepX on block 2'_`keepX on block 3')



}

\author{Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo.}

\seealso{ \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}} \code{\link{tune.block.splsda}} and http://www.mixOmics.org for more details.}

\examples{
\dontrun{
## validation for objects of class 'splsda'

data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)
out = tune.splsda(X, Y, ncomp = 3, nrepeat = 5, logratio = "none",
test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
progressBar = TRUE)


plot(out)
plot(out, sd=FALSE)


\dontrun{
## validation for objects of class 'mint.splsda'

data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

out = tune(method="mint.splsda", X=data,Y=type.id, ncomp=2, study=exp, test.keepX=seq(1,10,1))
out$choice.keepX

plot(out)



## validation for objects of class 'mint.splsda'

data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
protein = breast.TCGA$data.train$protein)
# set up a full design where every block is connected
# could also consider other weights, see our mixOmics manuscript
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = 5

# Tuning the first two components
# -------------

# definition of the keepX value to be tested for each block mRNA miRNA and protein
# names of test.keepX must match the names of 'data'
test.keepX = list(mrna = seq(10,40,20), mirna = seq(10,30,10), protein = seq(1,10,5))

# the following may take some time to run, note that for through tuning
# nrepeat should be > 1
tune = tune.block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
ncomp = ncomp, test.keepX = test.keepX, design = design, nrepeat = 3)

tune$choice.ncomp
tune$choice.keepX

plot(tune)
}
}
}

\keyword{regression}
\keyword{multivariate}
\keyword{hplot}
\value{none}
