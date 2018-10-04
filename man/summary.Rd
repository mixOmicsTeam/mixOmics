\name{summary}
\encoding{UTF-8}
\alias{summary}
\alias{summary.rcc}
\alias{summary.mixo_pls}
\alias{summary.mixo_spls}

\title{Summary Methods for CCA and PLS objects}

\description{
Produce \code{summary} methods for class \code{"rcc"},
\code{"pls"} and \code{"spls"}.
}

\usage{
\method{summary}{rcc}(object, what = c("all", "communalities", "redundancy"),
cutoff = NULL, digits = 4, \ldots)

\method{summary}{mixo_pls}(object, what = c("all", "communalities", "redundancy",
"VIP"), digits = 4, keep.var = FALSE, \ldots)

\method{summary}{mixo_spls}(object, what = c("all", "communalities", "redundancy",
"VIP"), digits = 4, keep.var = FALSE, \ldots)
}

\arguments{
\item{object}{object of class inherited from \code{"rcc"}, \code{"pls"} or \code{"spls"}.}
\item{cutoff}{real between 0 and 1. Variables with all correlations components below this cutoff
in absolute value are not showed (see Details).}
\item{digits}{integer, the number of significant digits to use when printing.
Defaults to \code{4}.}
\item{what}{character string or vector. Should be a subset of \code{c("all"},
\code{"summarised"}, \code{"communalities"}, \code{"redundancy"}, \code{"VIP"}).
\code{"VIP"} is only available for (s)PLS. See Details.}
\item{keep.var}{boolean. If \code{TRUE} only the variables with loadings not zero
(as selected by \code{spls}) are showed. Defaults to \code{FALSE}.}
\item{\ldots}{not used currently.}
}

\details{
The information in the \code{rcc}, \code{pls} or \code{spls} object is summarised, it includes:
the dimensions of \code{X} and \code{Y} data,
the number of variates considered, the canonical correlations (if \code{object} of class \code{"rcc"})
and the (s)PLS algorithm used (if \code{object} of class \code{"pls"} or \code{"spls"}) and
the number of variables selected on each of the sPLS components (if \code{x} of class \code{"spls"}).

\code{"communalities"} in \code{what} gives Communalities Analysis.
\code{"redundancy"} display Redundancy Analysis.
\code{"VIP"} gives the Variable Importance in the Projection (VIP)
coefficients fit by \code{pls} or \code{spls}.
If \code{what} is \code{"all"}, all are given.

For class \code{"rcc"}, when a value to \code{cutoff} is specified, the correlations
between each variable and the equiangular vector between
\eqn{X}- and \eqn{Y}-variates are computed. Variables with at least one correlation componente
bigger than \code{cutoff} are showed. The defaults is \code{cutoff=NULL} all the variables
are given.
}

\value{
The function \code{summary} returns a list with components:
\item{ncomp}{the number of components in the model.}
\item{cor}{the canonical correlations.}
\item{cutoff}{the cutoff used.}
\item{keep.var}{list containing the name of the variables selected.}
\item{mode}{the algoritm used in \code{pls} or \code{spls}.}
\item{Cm}{list containing the communalities.}
\item{Rd}{list containing the redundancy.}
\item{VIP}{matrix of VIP coefficients.}
\item{what}{subset of \code{c("all"}, \code{"communalities"}, \code{"redundancy"}, \code{"VIP"}).}
\item{digits}{the number of significant digits to use when printing.}
\item{method}{method used: \code{rcc}, \code{pls} or \code{spls}.}

}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{rcc}}, \code{\link{pls}}, \code{\link{spls}}, \code{\link{vip}}.}

\examples{
## summary for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
more <- summary(nutri.res, cutoff = 0.65)

\dontrun{
## summary for objects of class 'pls'
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)
more <- summary(linn.pls)

## summary for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))
more <- summary(toxicity.spls, what = "redundancy", keep.var = TRUE)
}
}

\keyword{regression}
\keyword{multivariate}
