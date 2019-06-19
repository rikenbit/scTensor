\name{convertToNCBIGeneID}
\docType{methods}
\alias{convertToNCBIGeneID}

\title{
ID conversion function to create the input matrix for scTensor
}
\description{
Duplicated row names or NA row names are removed from input matrix,
and the matrix with the row names of NCBI Gene ID is generated.
}
\usage{
convertToNCBIGeneID(input, rowID, LefttoRight)
}

\arguments{
\item{input}{
Input matrix (or data.frame).
}
\item{rowID}{
Gene identifier in each row of the input matrix.
The length of rowID must be same as the number of tht input matrix.
}
\item{LefttoRight}{
Corresponding table with left column (Gene identifier <same type of row ID>) and right column (NCBI Gene ID).
}
}

\value{
A matrix object with with the row names of NCBI Gene ID.
}

\author{Koki Tsuyuzaki}

\examples{
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", NA, "B", "B")
LefttoRight <- rbind(
  c("A", "1"),
  c("B", "2"),
  c("B", "4"),
  c("D", NA)
)
(input <- convertToNCBIGeneID(input, rowID, LefttoRight))
}

\keyword{methods}