\name{stemcells}
\encoding{UTF-8}
\alias{stemcells}

\docType{data}

\title{Human Stem Cells Data}

\description{
This data set contains the expression of a random subset of 400 genes in 125 samples from 4 independent studies and 3 cell types.
}

\usage{data(stemcells)}

\format{
A list containing the following components:
\describe{
\item{\code{gene}}{data matrix with 125 rows and 400 columns. Each row represents
an experimental sample, and each column a single gene.}
\item{\code{celltype}}{a factor indicating the cell type of each sample.}
\item{\code{study}}{a factor indicating the study from which the sample was extracted.}
}
}

\details{
This data set contains the expression of a random subset of 400 genes in 125 samples from 4 independent studies and 3 cell types. Those studies can be combined and analysed using the MINT procedure.
}


\references{
Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms. BMC Bioinformatics 18:128.
}

\keyword{datasets}
\value{none}
