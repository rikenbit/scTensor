\docType{methods}
\name{setParam}
\alias{setParam}
\alias{setParam<-}
\alias{setParam,CCSParams-method}
\alias{setParam<-,CCSParams-method}
\alias{setParam,CCSParams,ANY-method}
\alias{setParam<-,CCSParams,ANY-method}
\title{Set a parameter}
\usage{
setParam(object, name) <- value
\S4method{setParam}{CCSParams}(object, name, value)
}
\arguments{
\item{object}{object to set parameter in.}

\item{name}{name of the parameter to set.}

\item{value}{value to set the paramter to.}
}
\value{
Object with new parameter value.
}
\description{
Function for setting parameter values.
}
\examples{
params <- newCCSParams()

setParam(params, "nGene") <- 20000
setParam(params, "nCell") <- c(12, 43, 323)
setParam(params, "cciInfo") <- list(nPair=2000,
                                   CCI1=list(
                                       LPattern=c(1,0,0),
                                       RPattern=c(0,1,1),
                                       nGene=100,
                                       fc="E10"),
                                   CCI2=list(
                                       LPattern=c(0,0,1),
                                       RPattern=c(1,1,1),
                                       nGene=200,
                                       fc="E10"),
                                   CCI3=list(
                                       LPattern=c(1,1,1),
                                       RPattern=c(1,0,1),
                                       nGene=300,
                                       fc="E10")
                                   )
setParam(params, "lambda") <- 0.1
setParam(params, "seed") <- 111
}