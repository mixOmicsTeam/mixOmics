\name{plot.rcc}
\encoding{UTF-8}
\alias{plot.rcc}

\title{Canonical Correlations Plot}

\description{
This function provides scree plot of the canonical correlations.
}

\usage{
\method{plot}{rcc}(x, scree.type = c("pointplot", "barplot"), \ldots)
}

\arguments{
\item{x}{object of class inheriting from \code{"rcc"}.}
\item{scree.type}{character string, (partially) matching one of \code{"pointplot"} or \code{"barplot"},
determining the kind of	scree plots to be produced.}
\item{\ldots}{arguments to be passed to other methods. For the \code{"pointplot"} type
see \code{\link{points}}, for \code{"barplot"} type see \code{\link{barplot}}.}
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{points}}, \code{\link{barplot}}, \code{\link{par}}.}

\examples{
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)

## 'pointplot' type scree
plot(nutri.res) #(default)

\dontrun{
plot(nutri.res, pch = 19, cex = 1.2,
col = c(rep("red", 3), rep("darkblue", 18)))

## 'barplot' type scree
plot(nutri.res, scree.type = "barplot")

plot(nutri.res, scree.type = "barplot", density = 20, col = "black")
}
}

\keyword{multivariate}
\keyword{hplot}
\value{none}
