% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellData-class.R
\docType{class}
\name{cellData-class}
\alias{[,cellData,character,ANY-method}
\alias{cellData}
\alias{cellData,}
\alias{cellData-class}
\alias{cellData-method}
\alias{show,cellData-method}
\title{The cellData class}
\usage{
cellData(name, positive.ctr = character(0), negative.ctr = character(0),
  expwell = character(0),
  norm.method = getOption("opm.normalization.method"),
  QC.threshold = getOption("opm.QC.threshold"))

\S4method{[}{cellData,character,ANY}(x, i)

\S4method{show}{cellData}(object)
}
\arguments{
\item{name}{character, the analyzed item}

\item{positive.ctr}{a character vector, the positive control well IDs,
e.g. c("E05", "E06")}

\item{negative.ctr}{a character vector, the positive control well IDs,
e.g. c("B05", "B06")}

\item{expwell}{include all wells except control and neglect.well if NULL}

\item{norm.method}{character the normalization method.}

\item{QC.threshold}{numeric, the thresholds in the quality control.}

\item{x}{a cellData object}

\item{i}{a requested slot name}

\item{object}{a cellData class}
}
\value{
a cellData object
}
\description{
The main class used in \code{OperaMate} to hold all levels of experiment
data of a specific type.
}
\section{Slots}{

\describe{
\item{\code{name}}{character, one parameter in the Columbus system report.}

\item{\code{posctrwell}}{a character vector, the positive control well IDs, e.g. B05.}

\item{\code{negctrwell}}{a character vector, the negative control well IDs, e.g. B05.}

\item{\code{expwell}}{a character vector, the sample well IDs, e.g. C15.}

\item{\code{cellNum}}{matrix, cell numbers}

\item{\code{origin.data}}{a numeric matrix,
the raw data matrix with rows the well IDs and columns the plate IDs.}

\item{\code{norm.data}}{a numeric matrix, the normalized data.}

\item{\code{qc.data}}{a numeric matrix, the data after quality
control, with the rows are "barcode:wellID" and columns are the data of
all replicated samples and their means, and if they have passed the quality
control.}

\item{\code{norm.method}}{character the normalization method.}

\item{\code{QC.threshold}}{numeric, the thresholds in the quality control.}

\item{\code{plate.quality}}{a logical matrix, the quality data with
the rows are the barcode and columns are the replicateIDs.}

\item{\code{plate.quality.data}}{a list of plate correlations and plate z' factors}

\item{\code{Sig}}{a list of the following components:
\itemize{
\item{\code{SigMat}:}{a logic matrix marking the high and low expressed hits}
\item{\code{threshold}:}{the threshold of the high and low expressed hits}
\item{\code{stats}:}{the numbers of the high and low expressed hits}
\item{\code{pvalue}:}{the pvalue of each sample by t tests}
}}
}}
\section{Methods}{

\describe{
\item{Constructor}{
\code{cellData(name, positive.ctr = character(0), negative.ctr = character(0),
expwell = character(0), norm.method =  getOption("opm.normalization.method"),
QC.threshold = getOption("opm.QC.threshold")}.
}
\item{Show}{signature(object="cellLoad"). Displays object content as text.}
\item{Accessor}{x[i]. \code{x}: a cellData object;
\code{i}: character, a cellData slot name.}
}
}
\examples{
oneCell <- cellData(name = "Average Intensity of Nuclei",
positive.ctr = c("H02", "J02", "L02"),
negative.ctr = c("C23", "E23", "G23"))
oneCell
oneCell["name"]
}

