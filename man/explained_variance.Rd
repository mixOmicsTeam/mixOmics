\name{explained_variance}
\encoding{UTF-8}
\alias{explained_variance}

\title{Calculation of explained variance}

\description{
This function calculates the variance explained by variates.
}

\usage{
explained_variance(data, variates, ncomp)
}

\arguments{
\item{data}{numeric matrix of predictors}
\item{variates}{variates as obtained from a \code{pls} object for instance}
\item{ncomp}{number of components. Should be lower than the number of columns of \code{variates}}

}

\details{

\code{explained_variance} calculates the explained variance of each variates out of the total variance in \code{data}.
}

\value{
\code{explained_variance} simply returns the explained variance for each variate.
}


\author{Florian Rohart}

\seealso{\code{\link{spls}}, \code{\link{splsda}},
\code{\link{plotIndiv}}, \code{\link{plotVar}},
\code{\link{cim}}, \code{\link{network}}.}

\examples{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

ex = explained_variance(toxicity.spls$X, toxicity.spls$variates$X, ncomp =2)

# ex should be the same as
toxicity.spls$explained_variance$X

}

\keyword{regression}
\keyword{multivariate}
