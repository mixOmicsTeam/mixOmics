\name{study_split}
\encoding{UTF-8}
\alias{study_split}

\title{divides a data matrix in a list of matrices defined by a factor}

\description{
\code{study_spli} divides a data matrix in a list of matrices defined by a \code{study} input.
}

\usage{
study_split(data, study)
}

\arguments{
\item{data}{numeric matrix of predictors}
\item{study}{grouping factor indicating which samples are from the same study}

}

\value{
\code{study_split} simply returns a list of the same length as the number of levels of \code{study} tha contains submatrices of \code{data}.
}


\author{Florian Rohart}

\seealso{\code{\link{mint.pls}}, \code{\link{mint.spls}}, \code{\link{mint.plsda}},
\code{\link{mint.splsda}}.}

\examples{
data(stemcells)
data = stemcells$gene
exp = stemcells$study

data.list = study_split(data, exp)

names(data.list)
lapply(data.list, dim)
table(exp)
}


\keyword{regression}
\keyword{multivariate}
