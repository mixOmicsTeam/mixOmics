% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{network}
\alias{network}
\alias{network.default}
\alias{network.rcc}
\alias{network.pls}
\alias{network.spls}
\title{Relevance Network for (r)CCA and (s)PLS regression}
\usage{
network(
  mat,
  comp = NULL,
  blocks = c(1, 2),
  cutoff = 0,
  row.names = TRUE,
  col.names = TRUE,
  block.var.names = TRUE,
  graph.scale = 0.5,
  size.node = 0.5,
  color.node = NULL,
  shape.node = NULL,
  alpha.node = 0.85,
  cex.node.name = NULL,
  color.edge = color.GreenRed(100),
  lty.edge = "solid",
  lwd.edge = 1,
  show.edge.labels = FALSE,
  cex.edge.label = 1,
  show.color.key = TRUE,
  symkey = TRUE,
  keysize = c(1, 1),
  keysize.label = 1,
  breaks,
  interactive = FALSE,
  layout.fun = NULL,
  save = NULL,
  name.save = NULL,
  plot.graph = TRUE
)
}
\arguments{
\item{mat}{numeric matrix of values to be represented. Alternatively,
an object from one of the following models: \code{mix_pls}, \code{plsda}, 
\code{mixo_spls}, \code{splsda}, \code{rcc}, \code{sgcca}, \code{rgcca}, 
\code{sgccda}.}

\item{comp}{atomic or vector of positive integers. The components to
adequately account for the data association. Defaults to \code{comp = 1}.}

\item{blocks}{a vector indicating the block variables to display.}

\item{cutoff}{numeric value between \code{0} and \code{1}. The tuning
threshold for the relevant associations network (see Details).}

\item{row.names, col.names}{character vector containing the names of \eqn{X}-
and \eqn{Y}-variables.}

\item{block.var.names}{either a list of vector components for variable names
in each block or FALSE for no names. If TRUE, the columns names of the
blocks are used as names.}

\item{graph.scale}{Numeric between 0 and 1 which alters the scale of the entire
plot. Increasing the value decreases the size of nodes and increases their distance
from one another. Defaults to 0.5.}

\item{size.node}{Numeric between 0 and 1 which determines the relative size of nodes.
Defaults to 0.5.}

\item{color.node}{vector of length two, the colors of the \eqn{X} and
\eqn{Y} nodes (see Details).}

\item{shape.node}{character vector of length two, the shape of the \eqn{X}
and \eqn{Y} nodes (see Details).}

\item{alpha.node}{Numeric between 0 and 1 which determines the opacity of nodes.
Only used in block objects.}

\item{cex.node.name}{the font size for the node labels.}

\item{color.edge}{vector of colors or character string specifying the colors
function to using to color the edges, set to default to
\code{color.GreenRed(100)} but other palettes can be chosen (see Details and
Examples).}

\item{lty.edge}{character vector of length two, the line type for the edges
(see Details).}

\item{lwd.edge}{vector of length two, the line width of the edges (see
Details).}

\item{show.edge.labels}{logical. If \code{TRUE}, plot association values as
edge labels (defaults to \code{FALSE}).}

\item{cex.edge.label}{the font size for the edge labels.}

\item{show.color.key}{Logical. If \code{TRUE} a color key should be plotted.}

\item{symkey}{Logical indicating whether the color key should be made
symmetric about 0. Defaults to \code{TRUE}.}

\item{keysize}{numeric value indicating the size of the color key.}

\item{keysize.label}{vector of length 1, indicating the size of the labels
and title of the color key.}

\item{breaks}{(optional) either a numeric vector indicating the splitting
points for binning \code{mat} into colors, or a integer number of break
points to be used, in which case the break points will be spaced equally
between \code{min(mat)} and \code{max(mat)}.}

\item{interactive}{logical. If \code{TRUE}, a scrollbar is created to change
the cutoff value interactively (defaults to \code{FALSE}). See Details.}

\item{layout.fun}{a function. It specifies how the vertices will be placed
on the graph. See help(layout) in the igraph package. Defaults to
layout.fruchterman.reingold.}

\item{save}{should the plot be saved ? If so, argument to be set either to
\code{'jpeg'}, \code{'tiff'}, \code{'png'} or \code{'pdf'}.}

\item{name.save}{character string giving the name of the saved file.}

\item{plot.graph}{logical. If \code{TRUE} (default), plotting window will be 
filled with network. If \code{FALSE}, then no graph will be plotted, though 
the return value of the function is the exact same.}
}
\value{
\code{network} return a list containing the following components:
\item{M}{the correlation matrix used by \code{network}.} \item{gR}{a
\code{graph} object to save the graph for cytoscape use (requires to load
the \pkg{igraph} package).}
}
\description{
Display relevance associations network for (regularized) canonical
correlation analysis and (sparse) PLS regression. The function avoids the
intensive computation of Pearson correlation matrices on large data set by
calculating instead a pair-wise similarity matrix directly obtained from the
latent components of our integrative approaches (CCA, PLS, block.pls
methods). The similarity value between a pair of variables is obtained by
calculating the sum of the correlations between the original variables and
each of the latent components of the model. The values in the similarity
matrix can be seen as a robust approximation of the Pearson correlation (see
González et al. 2012 for a mathematical demonstration and exact formula).
The advantage of relevance networks is their ability to simultaneously
represent positive and negative correlations, which are missed by methods
based on Euclidean distances or mutual information. Those networks are
bipartite and thus only a link between two variables of different types can
be represented. The network can be saved in a .glm format using the
\code{igraph} package, the function \code{write.graph} and extracting the
output \code{object$gR}, see details.
}
\details{
\code{network} allows to infer large-scale association networks between the
\eqn{X} and \eqn{Y} datasets in \code{rcc} or \code{spls}. The output is a
graph where each \eqn{X}- and \eqn{Y}-variable corresponds to a node and the
edges included in the graph portray associations between them.

In \code{rcc}, to identify \eqn{X}-\eqn{Y} pairs showing relevant
associations, \code{network} calculate a similarity measure between \eqn{X}
and \eqn{Y} variables in a pair-wise manner: the scalar product value
between every pairs of vectors in dimension \code{length(comp)} representing
the variables \eqn{X} and \eqn{Y} on the axis defined by \eqn{Z_i} with
\eqn{i} in \code{comp}, where \eqn{Z_i} is the equiangular vector between
the \eqn{i}-th \eqn{X} and \eqn{Y} canonical variate.

In \code{spls}, if \code{object$mode} is \code{regression}, the similarity
measure between \eqn{X} and \eqn{Y} variables is given by the scalar product
value between every pairs of vectors in dimension \code{length(comp)}
representing the variables \eqn{X} and \eqn{Y} on the axis defined by
\eqn{U_i} with \eqn{i} in \code{comp}, where \eqn{U_i} is the \eqn{i}-th
\eqn{X} variate. If \code{object$mode} is \code{canonical} then \eqn{X} and
\eqn{Y} are represented on the axis defined by \eqn{U_i} and \eqn{V_i}
respectively.

Variable pairs with a high similarity measure (in absolute value) are
considered as relevant. By changing the cut-off, one can tune the relevance
of the associations to include or exclude relationships in the network.

\code{interactive=TRUE} open two device, one for association network, one
for scrollbar, and define an interactive process: by clicking either at each
end (\eqn{-} or \eqn{+}) of the scrollbar or at middle portion of this.
The position of the slider indicate which is the `cutoff' value associated
to the display network.

The network can be saved in a .glm format using the \pkg{igraph} package,
the function \code{write.graph} and extracting the output \code{obkect$gR}.

The interactive process is terminated by clicking the second button and
selecting \code{Stop} from the menu, or from the \code{Stop} menu on the graphics
window.

The \code{color.node} is a vector of length two, of any of the three kind of
\code{R} colors, i.e., either a color name (an element of \code{colors()}),
a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i}
meaning \code{palette()[i]}. \code{color.node[1]} and \code{color.node[2]}
give the color for filled nodes of the \eqn{X}- and \eqn{Y}-variables
respectively. Defaults to \code{c("white", "white")}.

\code{color.edge} give the color to edges with colors corresponding to the
values in \code{mat}. Defaults to \code{color.GreenRed(100)} for negative
(green) and positive (red) correlations. We also propose other palettes of
colors, such as \code{color.jet} and \code{color.spectral}, see help on
those functions, and examples below. Other palette of colors from the stats
package can be used too.

\code{shape.node[1]} and \code{shape.node[2]} provide the shape of the nodes
associate to \eqn{X}- and \eqn{Y}-variables respectively. Current acceptable
values are \code{"circle"} and \code{"rectangle"}. Defaults to
\code{c("circle", "rectangle")}.

\code{lty.edge[1]} and \code{lty.egde[2]} give the line type to edges with
positive and negative weight respectively. Can be one of \code{"solid"},
\code{"dashed"}, \code{"dotted"}, \code{"dotdash"}, \code{"longdash"} and
\code{"twodash"}. Defaults to \code{c("solid", "solid")}.

\code{lwd.edge[1]} and \code{lwd.edge[2]} provide the line width to edges
with positive and negative weight respectively. This attribute is of type
double with a default of \code{c(1, 1)}.
}
\section{Warning}{
 If the number of variables is high, the generation of the
network generation can take some time.
}

\examples{
## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

\dontrun{
# may not work on the Linux version, use Windows instead
# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example1-network.jpeg', res = 600, width = 4000, height = 4000)
network(nutri.res, comp = 1:3, cutoff = 0.6)
dev.off()

## Changing the attributes of the network

# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example2-network.jpeg')
network(nutri.res, comp = 1:3, cutoff = 0.45,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("circle", "rectangle"),
color.edge = color.jet(100),
lty.edge = "solid", lwd.edge = 2,
show.edge.labels = FALSE)
dev.off()


## interactive 'cutoff' - select the 'cutoff' and "see" the new network
## only run this during an interactive session
if (interactive()) {
    network(nutri.res, comp = 1:3, cutoff = 0.55, interactive = TRUE)
}
dev.off()

## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

# sometimes with Rstudio might not work because of margin issues,
# in that case save it as an image
jpeg('example3-network.jpeg')
network(toxicity.spls, comp = 1:3, cutoff = 0.8,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("rectangle", "circle"),
color.edge = color.spectral(100),
lty.edge = "solid", lwd.edge =  1,
show.edge.labels = FALSE, interactive = FALSE)
dev.off()
}
}
\references{
Mathematical definition: González I., Lê Cao K-A., Davis, M.J.
and Déjean, S. (2012). Visualising associations between paired omics data
sets. J. Data Mining 5:19.
\url{http://www.biodatamining.org/content/5/1/19/abstract}

Examples and illustrations:

Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
feature selection and multiple data integration. PLoS Comput Biol 13(11):
e1005752

Relevance networks:

Butte, A. J., Tamayo, P., Slonim, D., Golub, T. R. and Kohane, I. S. (2000).
Discovering functional relationships between RNA expression and
chemotherapeutic susceptibility using relevance networks. \emph{Proceedings
of the National Academy of Sciences of the USA} \bold{97}, 12182-12186.

Moriyama, M., Hoshida, Y., Otsuka, M., Nishimura, S., Kato, N., Goto, T.,
Taniguchi, H., Shiratori, Y., Seki, N. and Omata, M. (2003). Relevance
Network between Chemosensitivity and Transcriptome in Human Hepatoma Cells.
\emph{Molecular Cancer Therapeutics} \bold{2}, 199-205.
}
\seealso{
\code{\link{plotVar}}, \code{\link{cim}},
\code{\link{color.GreenRed}}, \code{\link{color.jet}},
\code{\link{color.spectral}} and http: //www.mixOmics.org for more details.
}
\author{
Ignacio González, Kim-Anh Lê Cao, AL J Abadi
}
\keyword{dplot}
\keyword{graphs}
\keyword{hplot}
\keyword{iplot}
\keyword{multivariate}
