% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellData-class.R
\docType{methods}
\name{cellNumLoad}
\alias{cellNumLoad}
\alias{cellNumLoad,cellData,cellData-method}
\title{Load cell number}
\usage{
cellNumLoad(object, object.cellnum)

\S4method{cellNumLoad}{cellData,cellData}(object, object.cellnum)
}
\arguments{
\item{object}{a cellData object}

\item{object.cellnum}{a cellData object for cell numbers}
}
\value{
a \code{cellData} object, with initialized slot \code{cellNum}
}
\description{
Load cell number
}
\examples{
data(demoCell)
data(demoCellNum)
oneCell <- cellNumLoad(oneCell, oneCellNum)
str(oneCell["cellNum"])
}

