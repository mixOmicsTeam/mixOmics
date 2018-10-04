\name{logratio.transfo}
\encoding{UTF-8}
\alias{logratio.transfo}

\title{Log-ratio transformation}

\description{
This function applies a log transformation to the data, either CLR or ILR
}

\usage{
logratio.transfo(X, logratio = "none", offset = 0)
}

\arguments{
\item{X}{numeric matrix of predictors}
\item{logratio}{log-ratio transform to apply, one of "none", "CLR" or "ILR"}
\item{offset}{Value that is added to X for CLR and ILR log transformation. Default to 0.}

}

\details{

\code{logratio.transfo} applies a log transformation to the data, either CLR (centered log ratio transformation) or ILR (Isometric Log Ratio transformation). In the case of CLR log-transformation, X needs to be a matrix of non-negative values and \code{offset} is used to shift the values away from 0, as commonly done with counts data.
}



\references{
Kim-Anh Lê Cao, Mary-Ellen Costello, Vanessa Anne Lakis, Francois Bartolo, Xin-Yi Chua, Remi Brazeilles, Pascale Rondeau
mixMC: a multivariate statistical framework to gain insight into Microbial Communities
bioRxiv 044206; doi: http://dx.doi.org/10.1101/044206

John Aitchison. The statistical analysis of compositional data. Journal of the Royal Statistical Society.
Series B (Methodological), pages 139-177, 1982.

Peter Filzmoser, Karel Hron, and Clemens Reimann. Principal component analysis for
compositional data with outliers. Environmetrics, 20(6):621-632, 2009.
}

\value{
\code{logratio.transfo} simply returns the log-ratio transformed data.
}


\author{Florian Rohart}

\seealso{\code{\link{pca}}, \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}}.}

\examples{
data(diverse.16S)
CLR = logratio.transfo(X = diverse.16S$data.TSS, logratio = 'CLR')
# no offset needed here as we have put it prior to the TSS, see www.mixOmics.org/mixMC

}
