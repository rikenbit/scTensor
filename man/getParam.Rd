\docType{methods}
\name{getParam}
\alias{getParam}
\alias{getParam,CCSParams-method}
\title{Get a parameter}
\usage{
getParam(object, name)

\S4method{getParam}{CCSParams}(object, name)
}
\arguments{
\item{object}{object to get parameter from.}

\item{name}{name of the parameter to get.}
}
\value{
The extracted parameter value
}
\description{
Accessor function for getting parameter values.
}
\examples{
params <- newCCSParams()

getParam(params, "nGene")
getParam(params, "nCell")
getParam(params, "cciInfo")
getParam(params, "lambda")
getParam(params, "seed")
}