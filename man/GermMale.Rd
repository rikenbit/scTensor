\name{GermMale}
\alias{GermMale}
\docType{data}
\title{
The matrix which is used as test data of scTensor.
}
\description{
A matrix with 1265 rows (genes) * 852 columns (cells).
}
\usage{data(GermMale)}
\details{
The data matrix is downloaded from GEO Series GSE86146
(https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86146&format=file).
Only male data is extracted and then the gene symbol is converted to
NCBI Gene ID by Homo.sapiens package.

For saving the package size, the number of genes are strictlly reduced by the
standard of highlly variable genes with threshold of p-value is 1E-50.
}
\references{
Li L. and Dong J. and Yan L. and Yong J. et al. (2017)
Single-Cell RNA-Seq Analysis Maps Development of Human Germline Cells
and Gonadal Niche Interactions. \emph{Cell Stem Cell}, \bold{20(6)}: 858-873
}
\seealso{
\code{\link{labelGermMale}}, \code{\link{tsneGermMale}}.
}
\examples{
data(GermMale)
}
\keyword{datasets}