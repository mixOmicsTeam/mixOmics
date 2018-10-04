\name{mat.rank}
\encoding{UTF-8}
\alias{mat.rank}

\title{Matrix Rank}

\description{
This function estimate the rank of a matrix.
}

\usage{
mat.rank(mat, tol)
}

\arguments{
\item{mat}{a numeric matrix or data frame that can contain missing values.}
\item{tol}{positive real, the tolerance for singular values, only those
with values larger than \code{tol} are considered non-zero.}
}

\details{
\code{mat.rank} estimate the rank of a matrix by
computing its singular values \eqn{d[i]} (using \code{nipals}).
The rank of the matrix can be defined as
the number of singular values \eqn{d[i] > 0}.

If \code{tol} is missing, it is given by
\code{tol=max(dim(mat))*max(d)*.Machine$double.eps}.
}

\value{
The returned value is a list with components:
\item{rank}{a integer value, the matrix rank.}
\item{tol}{the tolerance used for singular values.}
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{nipals}}}

\examples{
## Hilbert matrix
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
mat <- hilbert(16)
mat.rank(mat)

\dontrun{
## Hilbert matrix with missing data
idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
m.na <- m <- hilbert(9)[, 1:6]
m.na[idx.na == 0] <- NA
mat.rank(m)
mat.rank(m.na)
}
}


\keyword{algebra}


