#' Data extraction from one report
#'
#' Extracts data in the report to the slot \code{data} in the \code{expData}
#' object. An inner function of \code{\link{loadAll}}.
#'
#' @param onePlate an expData object
#' @param well.digits the digits of the well column in the well-gene
#' specification file
#' 
#' @return an \code{expData} object with initialized slot \code{data}.
#' @docType methods
#' @examples
#' datapath <- file.path(system.file("Test", package = "OperaMate"), "Tab")
#' lstPlates <- loadAll(cellformat = "Tab", datapath = datapath )
#' onePlate <- parseTemplete(lstPlates[[1]])
#' @export
parseTemplete <- function(onePlate, well.digits = 2) {
    if (onePlate["format"] == "Tab") {
        x <- parseTab(onePlate, well.digits)
    } else if (onePlate["format"] == "Matrix") {
        x <- parseMatrix(onePlate, well.digits)
    } else if (onePlate["format"] == "Full-Table") {
        x <- parseFullTab(onePlate, well.digits)
    } else {
        warning("Format is not supported. Cannot load the object Data.")
    }
    return(x)
}

## Data extraction when the report is of table format
##
## Extracts data in the report which is of table format. An inner function
## of \code{\link{parseTemplete}}.
##
## @param onePlate an expData object
## @return an \code{expData} object with initialized slot \code{data}.
## @docType methods
parseTab <- function(onePlate, well.digits = 2) {
    if (onePlate["format"] != "Tab") {
        stop("Error processor for ",
             onePlate["name"], ": ", onePlate["format"], ".")
    }
    expdata <- read.delim(onePlate["path"],
                          sep = "\t", stringsAsFactors = FALSE, skip = 3)
    if (colnames(expdata)[1] != "Row") {
        stop("Error processor for ",
             onePlate["name"], ": ", onePlate["format"], ".")
    }
    wellname <- paste0(expdata[, "Row"],
                       formatC(as.numeric(expdata[, "Column"]),
                               width = well.digits, flag = "0"))
    twellID <- ( !is.na(expdata[, 3]) & expdata[, 3] != "")
    expdata <- expdata[twellID, ]
    dataLoad(onePlate, as.list(expdata[, 4:ncol(expdata)]), wellname[twellID])
}

## Data extraction when the report is of table-full format
##
## Extracts data in the report which is of table format. An inner function
## of \code{\link{parseTemplete}}.
##
## @param onePlate an expData object
## @return an \code{expData} object with initialized slot \code{data}.
## @docType methods
parseFullTab <- function(onePlate, well.digits = 2) {
    if (onePlate["format"] != "Full-Table") {
        stop("Error processor for ",
             onePlate["name"], ": ", onePlate["format"], ".")
    }
    expdata <- read.delim(onePlate["path"],
                          sep = "\t", stringsAsFactors = FALSE)
    wellname <- paste0(LETTERS[expdata[, "Row"]],
                       formatC(as.numeric(expdata[, "Column"]),
                               width = well.digits, flag = "0"))
    dataLoad(onePlate, expdata, wellname)
}


## Data extraction when the report is of matrix format
##
## Extracts data in the report which is of matrix format. An inner function
## of \code{parseTemplete}.
##
## @param onePlate an expData object
## @return an \code{expData} object with initialized slot \code{data}.
## @docType methods
parseMatrix <- function(onePlate, well.digits = 2) {
    if (onePlate["format"] != "Matrix") {
        stop("Error processor for ",
             onePlate["name"], ": ", onePlate["format"], ".")
    }
    expdata <- read.delim(onePlate["path"],
                          sep = "\t",stringsAsFactors = FALSE, skip = 3,
                          header = FALSE)
    sline <- grep("Parameter: ", expdata[, 1])
    ## Format Checking
    region <- expdata[(sline[1] + 2):(sline[2] - 2), 2:ncol(expdata)]
    rownames(region) <- expdata[(sline[1] + 2):(sline[2] - 2), 1]
    colnames(region) <- formatC(seq(1, ncol(region)),
                                width = well.digits, flag = "0")
    vecregion <- as.vector(as.matrix(region))
    wellname <- apply( expand.grid(rownames(region), colnames(region))
                      , 1, paste, collapse = "")
    twellID <- ( !is.na(as.vector(as.matrix(region)))
       & as.vector(as.matrix(region)) != "")
    lstmatrix <- lapply(sline[2:length(sline)], function(line) {
        fmatrix <-
          expdata[(line + 2):(line + 1 + nrow(region)), 2:ncol(expdata)]
        f.3m <- as.vector(as.matrix(fmatrix))
        mode(f.3m) <- "numeric"
        return(f.3m[twellID])
    })
    wellname <- wellname[twellID]
    names(lstmatrix) <- gsub(" ", ".", expdata[sline[2:length(sline)], 2])
    dataLoad(onePlate, lstmatrix, wellname)
}   

#' Data importing
#' 
#' Initializes a list of \code{expData} objects from the Columbus system
#' reports.
#'
#' To facility the automatic file name parsing, the reports obtained from
#' Columbus system should be of the same format, and located under the same
#' directory. Users can obtain this plate specification table for further
#' modification. An example of the table can be referred by \code{\link{platemap}}.
#' After modification, users can submit a plate speficication data frame
#' to parameter \code{platemap}.
#' The data format supported for the reports are "Tab" and "Matrix".
#' If the reports are of other cellformats, you can specify
#' its cellformat and rewrite the function \code{parseTemplete} to import the
#' data seperately.
#' 
#' @param cellformat character specifying the format of the reports. Enable when
#' \code{platemap} is NULL.
#' @param datapath character specifying the location of the reports. Enable when
#' \code{platemap} is NULL.
#' @param egFilename a file name example
#' @param platemap data frame. See an example as \code{\link{platemap}}.
#' @param well.digits the digits of the well column in the well-gene
#' @return a list of \code{expData} objects
#'
#' @details An example of egFilename = list(eg.filename = "0205-s2-01.txt",
#' rep.id = "s2", exp.id = "01", sep = "-", barcode = "DSIMGA01").
#' well.digits: In the well-gene specification file, if the well ID is
#' B1, B2, ..., B11, the well.digit = 1; while B01, B02, ..., B11,
#' the well.digit = 2; and B001, B002, ..., B011, the well.digit =3.
#' @docType methods
#' @examples
#' # Data frame \code{platemap} provided
#' data(platemap)
#' platemap$Path <- file.path(
#' system.file("Test", package = "OperaMate"), platemap$Path)
#' lstPlates <- loadAll(platemap = platemap)
#' #
#' # Consistent file name format
#' datapath <- file.path(system.file("Test", package = "OperaMate"), "Tab")
#' egFilename <- list(eg.filename = "Tab.130504-s1-01.txt",
#' rep.id = "s1", exp.id = "01", sep = "-",
#' barcode = "DSIMGA01")
#' lstPlates <- loadAll(cellformat = "Tab", datapath = datapath,
#' egFilename = egFilename, well.digits = 2)
#' #
#' lstPlates[[1]]
#' @export
#'
loadAll <- function(cellformat = NULL, datapath = "./",
                    egFilename = getOption("opm.filename.example"),
                    well.digits = 2, platemap = NULL) {
    if (is.null(cellformat)) {
        ## Check if the existence and legality of basic information
        if (is.null(platemap)) {
            stop("Please specify data format or plate configuration table.
                 Check an example from data(platemap).")
        }
        if (any(! c("Barcode", "RepID" , "Format", "Path")
                %in% colnames(platemap) )) {
            stop("The format of mapping file is illegal! 
Please refer to data(platemap.example) for an example.")
        }
        if (any(! platemap$Format %in% c("Tab", "Matrix") )) {
            warning("Please choose a proper format from (Tab, Matrix).
Otherwise, you need to redefine the function parseTemplete to import data.")
        }
    } else {
        ## File names formats are consistent
        if (!cellformat %in% c("Tab", "Matrix")) {
            warning("Please choose a proper format from (Tab, Matrix).
Otherwise, you need to redefine the function parseTemplete to import data.")
        }
        if (is.null(egFilename)) {
            egFilename <- getOption("opm.filename.example")
        } else {
            tmp <- getOption("opm.filename.example")
            tmpID <- match(names(egFilename), names(tmp))
            tmp[!is.na(tmpID)] <- egFilename[na.omit(tmpID)]
            egFilename <- tmp
        }
        pfiles <- list.files(datapath, recursive = TRUE)
        if (length(pfiles) == 0)
          stop("No files in the provided file location.")
        platemap <- nameParser(pfiles, egFilename)
        platemap$Path <- file.path(datapath, pfiles)
        platemap$Format <- cellformat
    }

    lstPlates <- apply(platemap, 1, function(x) {
        onePlate <-
          expData( name  = paste(x["Barcode"], x["RepID"], sep = "-"),
                  path   = x["Path"],
                  rep.id = x["RepID"],
                  exp.id = x["Barcode"],
                  format = x["Format"])
                return(parseTemplete(onePlate, well.digits))
    })
    return(lstPlates)
}

