## Classes: expData:  data of each plate;

#' The expData class
#'
#' The expData class is a container to store data imported from one
#' Columbus system report
#' @rdname expData-class
#' @aliases expData
#' @aliases expData-class
#' @slot name character, the plate ID (barcode-replicateID), e.g. DSIMGA03-s1.
#' @slot path character, the path of the Columbus system report.
#' @slot rep.id character, replicateID, e.g. s1.
#' @slot exp.id character, barcode, e.g. DSIMGA03.
#' @slot data a list of vectors,
#' the vectorized raw data matrix of one plate of each type.
#' @slot format character, format of the Columbus system report.
#' @slot wellID a character vector, the well IDs.
#' @section Methods: 
#' \describe{
#' \item{Constructor}{
#' \code{expData(name, path, rep.id, exp.id, format)}.
#' }
#' \item{Show}{signature(object = "expData"). Displays object content as text.}
#' \item{Accessor}{x[i]. \code{x}: an expData object;
#' \code{i}: character, an expData slot name.}
#' \item{dataLoad}{\code{dataLoad(object, data, wellID)}}
#' }
#' @examples
#' onePlate <- expData(name = "130504-s1-02.txt",
#'          path = file.path(system.file("Test", package = "OperaMate"),
#'                            "Matrix", "130504-s1-02.txt"),
#'        rep.id = "s1",
#'        exp.id = "DSIMGA02",
#'        format = "Matrix")
#' onePlate
#' onePlate["name"]
#' @exportClass expData
#' @name expData-class
#' @aliases expData, expData-class
#' @rdname expData-class
setClass("expData", slots = c(name = "character",
                     path = "character",
                   rep.id = "character",
                   exp.id = "character",
                     data = "list",
                   format = "character",
                   wellID = "character"),
         prototype = list(format = c("Tab", "Matrix"))
         )

#' Constructor method of expData class.
#' @name expData
#' @rdname expData-class
#' @aliases expData, expData-method
#' @param name character, the plate ID (barcode-replicateID), e.g. DSIMGA03-s1.
#' @param path character, the path of the Columbus system report.
#' @param rep.id character, replicateID, e.g. s1.
#' @param exp.id character, barcode, e.g. DSIMGA03.
#' @param format character, format of the Columbus system report.
#' @return an expData object
#' @export
expData <- function(name, path, rep.id, exp.id, format) {
    if (!file.exists(path)) {
        stop("Cannot find file path:", path, ".")
    }
    if (length(format)>1) {
        format <- format[1]
    }

    new("expData", name = name, path = path,
        rep.id = rep.id, exp.id = exp.id, format = format)
}

#' Show method
#' @rdname expData-class
#' @param object a expData class
#' @aliases show,expData-method
setMethod("show", signature = (object = "expData"), function(object) {
    cat("An object of expData class.\n")
    cat("Name: ", object@name, ".\n", sep="")
    cat("Path: ", object@path, ".\n", sep="")
    cat("Data: ")
    str(object@data)
    cat("Wells: ")
    str(object@wellID)
})

## Getter
#' @rdname expData-class
#' @aliases [, expData-method
#' @param x a expData object
#' @param i a requested slot name
setMethod("[", signature = "expData", function(x, i) {
    if (i == "name")   { return(x@name)    } else {}
    if (i == "path")   { return(x@path)    } else {}
    if (i == "rep.id") { return(x@rep.id)  } else {}
    if (i == "exp.id") { return(x@exp.id)  } else {}
    if (i == "format") { return(x@format)  } else {}
    if (i == "wellID") { return(x@wellID)  } else {}
    if (i == "data")   { return(x@data)    } else {}
})

#' @name dataLoad
#' @rdname expData-class
#' @aliases dataLoad, expData-method
#' @param data the vectorized raw data matrix of one plate of each type.
#' @param wellID a character vector, the well IDs.
#' @export
setGeneric("dataLoad", function(object, data, wellID) {
    standardGeneric("dataLoad")
})
#' @rdname expData-class
#' @aliases dataLoad, expData-method
#' 
setMethod("dataLoad", signature = "expData", function(object, data, wellID) {
    if (!is.list(data)) {
        stop("'data' is a list with components the values of each parameter.")
    } else {}
    if (any( sapply(data, length) != length(wellID) )) {
        stop("Lengths of 'data' and 'wellID' do not match.")
    } else {}
    object@data <- data
    object@wellID <- wellID
    return(object)
})
