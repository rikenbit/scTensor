\name{tsneGermMale}
\alias{tsneGermMale}
\docType{data}
\title{
The result of Rtsne against GermMale
}
\description{
A List contains some parameters and the result of Rtsne function.
}
\usage{data(tsneGermMale)}
\details{
Rtsne is performed as follows.

library(Rtsne)
set.seed(123)
tsneGermMale <- Rtsne(dist(t(GermMale)), is_distance=TRUE, perplexity=40)
}
\references{
Li L. and Dong J. and Yan L. and Yong J. et al. (2017) Single-Cell RNA-Seq
Analysis Maps Development of Human Germline Cells and Gonadal Niche
Interactions. \emph{Cell Stem Cell}, \bold{20(6)}: 858-873
}
\seealso{
    \code{\link{labelGermMale}}, \code{\link{GermMale}}.
}
\examples{
data(tsneGermMale)
}
\keyword{datasets}