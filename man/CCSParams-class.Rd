\name{CCSParams-class}
\docType{class}
\alias{CCSParams-class}
\title{Class "CCSParams"}

\description{
The parameter object to be specified against cellCellSimulate function.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CCSParams", ...)}.
}

\section{Slots}{
   \describe{
    \item{nGene:}{The number of genes.}
    \item{nCell:}{The number of cells.}
    \item{cciInfo:}{The parameter to describe the CCI.}
    \item{lambda:}{The parameter for dropout simulation.}
    \item{seed:}{The seed for using random numbers.}
  }
}

\section{Methods}{
  \describe{
    \item{newCCSParams}{Generator of CCSParams object.}
    \item{getParam}{Getter function of the slot in CCSParams object.}
    \item{setParam<-}{Setter function of the slot in CCSParams object.}
  }
}

\seealso{
  \code{\link{newCCSParams}}, \code{\link{getParam}}, \code{\link{setParam<-}}
}

\keyword{classes}
