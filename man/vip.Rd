\name{vip}
\encoding{UTF-8}
\alias{vip}

\title{Variable Importance in the Projection (VIP)}

\description{
The function \code{vip} computes the influence on the \eqn{Y}-responses
of every predictor \eqn{X} in the model.
}

\usage{
vip(object)
}

\arguments{
\item{object}{object of class inheriting from \code{"pls"}, \code{"plsda"}, \code{"spls"} or  \code{"splsda"}.}
}

\details{
Variable importance in projection (VIP) coefficients reflect the relative
importance of each \eqn{X} variable for each \eqn{X} variate in the prediction
model. VIP coefficients thus represent
the importance of each \eqn{X} variable in fitting both the \eqn{X}-
and \eqn{Y}-variates,
since the \eqn{Y}-variates are predicted from the \eqn{X}-variates.

VIP allows to classify the \eqn{X}-variables according to their explanatory power of \eqn{Y}.
Predictors with large VIP, larger than 1, are the most relevant for explaining \eqn{Y}.
}

\value{
\code{vip} produces a matrix of VIP coefficients for each \eqn{X} variable (rows) on each
variate component (columns).
}

\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.
}

\author{Sébastien Déjean and Ignacio Gonz\`alez.}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{summary}}.}

\examples{
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)

linn.vip <- vip(linn.pls)

barplot(linn.vip,
beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan"),
ylim = c(0, 1.7), legend = rownames(linn.vip),
main = "Variable Importance in the Projection", font.main = 4)
}

\keyword{regression}
\keyword{multivariate}

