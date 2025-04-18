% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotArrow.R
\name{plotArrow}
\alias{plotArrow}
\title{Arrow sample plot}
\usage{
plotArrow(
  object,
  comp = c(1, 2),
  ind.names = TRUE,
  group = NULL,
  col = NULL,
  ind.names.position = c("start", "end"),
  ind.names.size = 2,
  pch = NULL,
  pch.size = 2,
  arrow.alpha = 0.6,
  arrow.size = 0.5,
  arrow.length = 0.2,
  legend = if (is.null(group)) FALSE else TRUE,
  legend.title = NULL,
  ...
)
}
\arguments{
\item{object}{object of class inheriting from \pkg{mixOmics}: \code{PLS,
sPLS, rCC, rGCCA, sGCCA, sGCCDA}}

\item{comp}{integer vector of length two (or three to 3d). The components
that will be used on the horizontal and the vertical axis respectively to
project the individuals.}

\item{ind.names}{either a character vector of names for the individuals to
be plotted, or \code{FALSE} for no names. If \code{TRUE}, the row names of
the first (or second) data matrix is used as names (see Details). If `pch` is set this
will overwrite the names as shapes. Default is \code{TRUE}. 
Not avaliable for MINT objects.}

\item{group}{Factor indicating the group membership for each sample.}

\item{col}{character (or symbol) color to be used. If \code{group} provided, should be a 
vector of the same length as groups, order of colours will be respected to correspond to 
order of \code{group} levels.}

\item{ind.names.position}{One of c('start', 'end') indicating where to show
the ind.names . Not used in block analyses, where centroids are used.}

\item{ind.names.size}{Numeric, sample name size.}

\item{pch}{plot character. A character string or a named vector of single
characters or integers whose names match those of \code{object$variates}.}

\item{pch.size}{Numeric, sample point character size.}

\item{arrow.alpha}{Numeric between 0 and 1 determining the opacity of arrows.}

\item{arrow.size}{Numeric, variable arrow head size.}

\item{arrow.length}{Numeric, length of the arrow head in 'cm'.}

\item{legend}{Logical, whether to show the legend if \code{group != NULL}.}

\item{legend.title}{Character, the legend title if \code{group != NULL}.}

\item{...}{Not currently used.
sample size to display sample names.}
}
\value{
A ggplot object
}
\description{
Represents samples from multiple coordinates to assess the alignment in the
latent space.
}
\details{
Graphical of the samples (individuals) is displayed in a superimposed manner
where each sample will be indicated using an arrow. The start of the arrow
indicates the location of the sample in \eqn{X} in one plot, and the tip the
location of the sample in \eqn{Y} in the other plot. Short arrows indicate a
strong agreement between the matching data sets, long arrows a disagreement
between the matching data sets. The representation space is scaled using the
range of coordinates so minimum and maximum values are equal for all blocks.
Since the algorithm maximises the covariance of these components, the
absolute values do not affect the alignment.

For objects of class \code{"GCCA"} and if there are more than 2 blocks, the
start of the arrow indicates the centroid between all data sets for a given
individual and the tips of the arrows the location of that individual in
each block.
}
\examples{

## 'spls' class - examples demonstrate how to control sample colours with sample names shown
# ------------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
spls.obj <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                 keepY = c(10, 10, 10))

# colors indicate time of necropsy, text is the dose, label at start of arrow
plotArrow(spls.obj,  group = as.factor(liver.toxicity$treatment[, 'Time.Group']),
          col = c("red", "blue", "purple", "darkgreen"), 
          ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
          legend = TRUE, position.names = 'start', legend.title = 'Time.Group')


## 'rcc' class - examples demonstrate how to control shape of all samples
# -------------------------------------------------------------------------------

# create model
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
rcc.obj <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

primary_groups <- nutrimouse$diet

# plot samples coloured by primary groups, by default shapes are all circles
plotArrow(rcc.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE)

# plot samples coloured by primary groups, force all samples to have 
# a different shape (2 = triangle)
plotArrow(rcc.obj, ind.names = FALSE,
          group = primary_groups, legend = TRUE,
          pch = 2)


## 'sgccda' class - examples demonstrate how to control shape of different blocks
# --------------------------------------------------------------------------------

data(breast.TCGA)
idx = seq(1, length(breast.TCGA$data.train$subtype), 10)
X <- list(mRNA = breast.TCGA$data.train$mrna[idx,],
          miRNA = breast.TCGA$data.train$mirna[idx,],
          protein = breast.TCGA$data.train$protein[idx,])
Y <- breast.TCGA$data.train$subtype[idx] # set the response variable

diablo.obj <- block.splsda(X, Y, ncomp = 2) # undergo multiblock sPLS-DA

# plot the samples using an arrow plot - 
plotArrow(diablo.obj, 
          ind.names = FALSE,
          legend = TRUE,
          title = 'TCGA, DIABLO comp 1 - 2') 

pchs <- c(3, 2, 1)
names(pchs) <- c("miRNA", "mRNA", "protein")

# plot the samples using an arrow plot changed shapes of data blocks
plotArrow(diablo.obj, 
          ind.names = FALSE,
          legend = TRUE,
          title = 'TCGA, DIABLO comp 1 - 2',
          pch = pchs) 
}
\references{
Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P.
(2009). Sparse canonical methods for biological data integration:
application to a cross-platform study. \emph{BMC Bioinformatics}
\bold{10}:34.
}
\seealso{
\code{\link{arrows}}, \code{\link{text}}, \code{\link{points}} and
http://mixOmics.org/graphics for more details.
}
\author{
Al J Abadi
}
\keyword{dplot}
\keyword{hplot}
\keyword{multivariate}
