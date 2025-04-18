% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLoadings.R, R/plotLoadings.mint.pls.R,
%   R/plotLoadings.mint.plsda.R, R/plotLoadings.pca.R, R/plotLoadings.pls.R,
%   R/plotLoadings.plsda.R
\name{plotLoadings}
\alias{plotLoadings}
\alias{plotLoadings.pls}
\alias{plotLoadings.spls}
\alias{plotLoadings.mint.pls}
\alias{plotLoadings.mint.spls}
\alias{plotLoadings.mint.plsda}
\alias{plotLoadings.mint.splsda}
\alias{plotLoadings.pca}
\alias{plotLoadings.mixo_pls}
\alias{plotLoadings.mixo_spls}
\alias{plotLoadings.rcc}
\alias{plotLoadings.sgcca}
\alias{plotLoadings.rgcca}
\alias{plotLoadings.mixo_plsda}
\alias{plotLoadings.mixo_splsda}
\alias{plotLoadings.sgccda}
\title{Plot of Loading vectors}
\usage{
plotLoadings(object, ...)

\method{plotLoadings}{mint.pls}(
  object,
  style = "graphics",
  comp = 1,
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.7,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  study = "global",
  ...
)

\method{plotLoadings}{mint.spls}(
  object,
  style = "graphics",
  comp = 1,
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.7,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  study = "global",
  ...
)

\method{plotLoadings}{mint.plsda}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.7,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  contrib = NULL,
  method = "mean",
  show.ties = TRUE,
  col.ties = "white",
  legend = TRUE,
  legend.color = NULL,
  legend.title = "Outcome",
  size.legend = 0.8,
  study = "global",
  ...
)

\method{plotLoadings}{mint.splsda}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.7,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  contrib = NULL,
  method = "mean",
  show.ties = TRUE,
  col.ties = "white",
  legend = TRUE,
  legend.color = NULL,
  legend.title = "Outcome",
  size.legend = 0.8,
  study = "global",
  ...
)

\method{plotLoadings}{pca}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  size.title = 2,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  ...
)

\method{plotLoadings}{mixo_pls}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = rel(1.5),
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  ...
)

\method{plotLoadings}{mixo_spls}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = rel(1.5),
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  ...
)

\method{plotLoadings}{rcc}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = rel(1.5),
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  ...
)

\method{plotLoadings}{sgcca}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = rel(1.5),
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  ...
)

\method{plotLoadings}{rgcca}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  col = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = rel(1.5),
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  ...
)

\method{plotLoadings}{mixo_plsda}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.6,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  contrib = NULL,
  method = "mean",
  show.ties = TRUE,
  col.ties = "white",
  legend = TRUE,
  legend.color = NULL,
  legend.title = "Outcome",
  size.legend = 0.8,
  ...
)

\method{plotLoadings}{mixo_splsda}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.6,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  contrib = NULL,
  method = "mean",
  show.ties = TRUE,
  col.ties = "white",
  legend = TRUE,
  legend.color = NULL,
  legend.title = "Outcome",
  size.legend = 0.8,
  ...
)

\method{plotLoadings}{sgccda}(
  object,
  comp = 1,
  style = "graphics",
  ndisplay = NULL,
  xlim = NULL,
  layout = NULL,
  border = NA,
  name.var = NULL,
  size.name = 0.7,
  title = NULL,
  subtitle,
  size.title = 2,
  size.subtitle = 1.6,
  size.axis = 0.7,
  X.label = NULL,
  Y.label = NULL,
  size.labs = 1,
  block,
  contrib = NULL,
  method = "mean",
  show.ties = TRUE,
  col.ties = "white",
  legend = TRUE,
  legend.color = NULL,
  legend.title = "Outcome",
  size.legend = 0.8,
  ...
)
}
\arguments{
\item{object}{object}

\item{\dots}{not used.}

\item{style}{argument to be set to either \code{'graphics'} or \code{'ggplot2'} 
to indicate the style of the plot. Default is \code{'graphics'}.}

\item{comp}{integer value, which component to plot. Default is 1.}

\item{ndisplay}{integer indicating how many of the most important variables
are to be plotted (ranked by decreasing weights in each component).
Useful to lighten a graph.}

\item{xlim}{Limit of the x-axis, vector of length 2.}

\item{layout}{Vector of two values (rows,cols) that indicates the layout of
the plot. If \code{layout} is provided, the remaining empty subplots are
still active. See details and examples.}

\item{col}{color of barplot, only for non-DA methods.}

\item{border}{Argument from \code{\link{barplot}}: indicates whether to draw
a border on the barplot and in which colour. Note that for -DA methods, 
when \code{contrib} is not provided, the colour is set to white and border black.}

\item{name.var}{A character vector or list indicating the names of the variables.
For single block analysis or when \code{study="all.partial"}, a vector of length equal 
to the number of variables in the block. For multi-block analysis, a list where each element 
is a vector of length equal to the number of variables in the corresponding block.
The names of the vector should match the names of the input data, see example.}

\item{size.name}{A numerical value giving the amount by which plotting the
variable name text should be magnified or reduced relative to the default.}

\item{title}{Title of the plot. Default is NULL.}

\item{subtitle}{subtitle for each plot, only used when several \code{block}
or \code{study} are plotted.}

\item{size.title}{size of the title}

\item{size.subtitle}{size of the subtitle}

\item{size.axis}{size of text on the X axis}

\item{X.label}{X axis label. Default is NULL.}

\item{Y.label}{Y axis label. Default is NULL.}

\item{size.labs}{size of the axis labels.}

\item{block}{A single value or vector indicating which block to plot. 
See details for behavior depending on object type.}

\item{study}{Indicates which study are to be plotted.  A character vector
containing some levels of \code{object$study}, "all.partial" to plot all
studies or "global" is expected.}

\item{contrib}{a character set to 'max' or 'min' indicating if the color of
the bar should correspond to the group with the maximal or minimal
expression levels / abundance.}

\item{method}{a character set to 'mean' or 'median' indicating the criterion
to assess the contribution. We recommend using median in the case of count
or skewed data.}

\item{show.ties}{Logical. If TRUE then tie groups appear in the color set by
\code{col.ties}, which will appear in the legend. Ties can happen when
dealing with count data type. By default set to TRUE.}

\item{col.ties}{Color corresponding to ties, only used if
\code{show.ties=TRUE} and ties are present.}

\item{legend}{Logical indicating if the legend indicating the group outcomes
should be added to the plot. Default value is TRUE.}

\item{legend.color}{A color vector of length the number of group outcomes.
See examples.}

\item{legend.title}{A set of characters to indicate the title of the legend.
Default value is NULL.}

\item{size.legend}{A numerical value giving the amount by which plotting the
legend text should be magnified or reduced relative to the default.}
}
\value{
Invisibly returns a \code{data.frame} containing the contribution of 
features on each component. For supervised models the contributions for
 each class is also specified. See details.
}
\description{
This function provides a horizontal bar plot to visualise loading vectors.
For discriminant analysis, it provides visualisation of highest or lowest
mean/median value of the variables with color code corresponding to the
outcome of interest.
}
\details{
The contribution of each variable for each component (depending on the
object) is represented in a barplot where each bar length corresponds to the
loading weight (importance) of the feature. The loading weight can be
positive or negative.

For discriminant analysis, the color corresponds to the group in which the
feature is most 'abundant'.  Note that this type of graphical output is
particularly insightful for count microbial data - in that latter case using
the \code{method = 'median'} is advised. Note also that if the parameter
\code{contrib} is not provided, plots are white.

For MINT analysis, \code{study="global"} plots the global loadings while
partial loadings are plotted when \code{study} is a level of
\code{object$study}. Since variable selection in MINT is performed at the
global level, only the selected variables are plotted for the partial
loadings even if the partial loadings are not sparse. See references.
Importantly for multi plots, the legend accounts for one subplot in the
layout design.

The \code{block} argument behavior varies depending on the object type. 
For \code{mixo_pls}, \code{mixo_spls}, \code{rcc}, \code{rgcca}, \code{sgcca},
\code{block} can be any number of blocks from \code{object$names$blocks}.
For code{mint.pls}, \code{mint.spls}: when \code{study="all.partial"}, 
can only be one block from \code{object$names$blocks}.
}
\examples{
## object of class 'spls'
# --------------------------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50),
keepY = c(10, 10))

plotLoadings(toxicity.spls)

# with xlim
xlim = c(-0.5, 0.5)
plotLoadings(toxicity.spls, xlim = xlim)


\dontrun{
## object of class 'splsda'
# --------------------------
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(paste0('treatment_' ,liver.toxicity$treatment[, 4]))

splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# contribution on comp 1, based on the median.
# Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 1, method = 'median')
plotLoadings(splsda.liver, comp = 1, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
#Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
# Colors indicate the group in which the median expression is minimal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# look at the contribution (median) for each variable
plot.contrib = plotLoadings(splsda.liver, comp = 2, method = 'median', plot = FALSE,
contrib = "max")
head(plot.contrib[[1]][,1:4])
# change the title of the legend and title name
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.title = 'Time',
title = 'Contribution plot', contrib = "max")

# no legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend = FALSE, contrib = "max")

# change the color of the legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4), contrib = "max")


# object 'splsda multilevel'
# -----------------

data(vac18)
X = vac18$genes
Y = vac18$stimulation
# sample indicates the repeated measurements
sample = vac18$sample
stimul = vac18$stimulation

# multilevel sPLS-DA model
res.1level = splsda(X, Y = stimul, ncomp = 3, multilevel = sample,
keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)

plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, size.name = 0.2, contrib = "max")

# too many transcripts? only output the top ones
plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, size.name = 0.5, ndisplay = 60, contrib = "max")



# object 'plsda'
# ----------------

# breast tumors
# ---
data(breast.tumors)
X = breast.tumors$gene.exp
Y = breast.tumors$sample$treatment

plsda.breast = plsda(X, Y, ncomp = 2)

name.var = as.character(breast.tumors$genes$name)
names(name.var) = colnames(X)

# with gene IDs, showing the top 60
plotLoadings(plsda.breast, contrib = 'max', comp = 1, method = 'median',
ndisplay = 60,
name.var = name.var,
size.name = 0.6,
legend.color = color.mixo(1:2))


# liver toxicity
# ---

data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$treatment[, 4]

plsda.liver = plsda(X, Y, ncomp = 2)

name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotLoadings(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100,
name.var = name.var, size.name = 0.4,
legend.color = color.mixo(1:4))


# object 'sgccda'
# ----------------

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda = wrapper.sgccda(X = data,
Y = Y,
design = design,
keepX = list(gene = c(10,10), lipid = c(15,15)),
ncomp = 2)

plotLoadings(nutrimouse.sgccda,block=2)
plotLoadings(nutrimouse.sgccda,block="gene")

}


#' PCA example
#' -----------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

#' Simple PCA plot
pca.liver = pca(X, ncomp = 2)
plotLoadings(pca.liver)

#' Customize PCA plot
plotLoadings(pca.liver, 
    comp = 2,
    ndisplay = 20,  # Show top 20 variables
    col = "steelblue",
    border = "black",
    size.name = 0.8,
    title = "PCA Loadings - Component 2",
    X.label = "Loading value",
    Y.label = "Variables")

#' PLS/SPLS example
#' ----------------
toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

#' Plot both blocks with custom layout
plotLoadings(toxicity.spls,
    comp = 2,
    block = c("X", "Y"),
    layout = c(2, 2),
    title = "PLS Loadings - Component 2",
    subtitle = c("Gene Block", "Clinical Block"))

#' PLSDA/SPLSDA example
#' -------------------
X = as.matrix(liver.toxicity$gene)
Y = as.factor(paste0('treatment_', liver.toxicity$treatment[, 4]))
splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))

#' Show contribution with gene names
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    name.var = name.var,
    size.name = 0.5,
    title = "Liver Treatment - Component 2",
    legend.title = "Treatment",
    legend.color = color.mixo(1:4))


#' MINT PLSDA example
#' -----------------
data(stemcells)
X = stemcells$gene
Y = stemcells$celltype
study = stemcells$study
mint.splsda = mint.splsda(X = X, Y = Y, ncomp = 3, keepX = c(10, 5, 15), study = study)

#' All partial loadings with custom layout
plotLoadings(mint.splsda,
    comp = 2,
    study = "all.partial",
    contrib = "max",
    method = "median",
    layout = c(2, 4),
    legend.color = color.mixo(1:3),
    subtitle = paste("Study", 1:4),
    show.ties = TRUE,
    col.ties = "gray")

#' Advanced customization examples
#' -----------------------------

#' Custom variable names and sizes
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    name.var = name.var,
    size.name = 0.8,
    title = "Custom Variable Names",
    size.title = 1.2,
    size.axis = 0.8,
    size.labs = 1.0)

#' No legend
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    legend = FALSE)

#' Custom legend
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    legend.title = "Treatment Groups",
    legend.color = c("red", "blue", "green", "purple"),
    size.legend = 0.8)

#' Minimal display with top variables
plotLoadings(splsda.liver,
    comp = 2,
    method = 'median',
    contrib = "max",
    ndisplay = 20,  # Show only top 20 variables
    title = "Top 20 Contributing Variables")

}
\references{
Rohart F. et al (2016, submitted). MINT: A multivariate
integrative approach to identify a reproducible biomarker signature across
multiple experiments and platforms.

Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group
PLS Regression: Application to Epidemiology. In New Perspectives in Partial
Least Squares and Related Methods, pages 243-255. Springer.

Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
molecular drivers from multi-omics assays, Bioinformatics, 
Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.

Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
Sparse canonical methods for biological data integration: application to a
cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.

Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
Editions Technic.

Wold H. (1966). Estimation of principal components and related models by
iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
Analysis}. Academic Press, N.Y., 391-420.
}
\seealso{
\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
\code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
\code{\link{mint.plsda}}, \code{\link{mint.splsda}},
\code{\link{block.pls}}, \code{\link{block.spls}},
\code{\link{block.plsda}}, \code{\link{block.splsda}},
\code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
\code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}
}
\author{
Florian Rohart, Kim-Anh Lê Cao, Benoit Gautier, Al J Abadi
}
\keyword{multivariate}
