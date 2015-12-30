## cellData: data of each type.
##
#' The cellData class
#'
#' The main class used in \code{OperaMate} to hold all levels of experiment
#' data of a specific type.
#' 
#' @slot name character, one parameter in the Columbus system report.
#' @slot posctrwell a character vector, the positive control well IDs, e.g. B05.
#' @slot negctrwell a character vector, the negative control well IDs, e.g. B05.
#' @slot expwell a character vector, the sample well IDs, e.g. C15.
#' @slot cellNum matrix, cell numbers
#' @slot origin.data a numeric matrix,
#' the raw data matrix with rows the well IDs and columns the plate IDs.
#' @slot norm.data a numeric matrix, the normalized data.
#' @slot qc.data a numeric matrix, the data after quality
#' control, with the rows are "barcode:wellID" and columns are the data of
#' all replicated samples and their means, and if they have passed the quality
#' control.
#' @slot norm.method character the normalization method.
#' @slot QC.threshold numeric, the thresholds in the quality control.
#' @slot plate.quality a logical matrix, the quality data with
#' the rows are the barcode and columns are the replicateIDs.
#' @slot plate.quality.data a list of plate correlations and plate z' factors
#' @slot Sig a list of the following components:
#' \itemize{
#' \item{\code{SigMat}:}{a logic matrix marking the high and low expressed hits}
#' \item{\code{threshold}:}{the threshold of the high and low expressed hits}
#' \item{\code{stats}:}{the numbers of the high and low expressed hits}
#' \item{\code{pvalue}:}{the pvalue of each sample by t tests}
#' }
#' 
#' @section Methods:
#' \describe{
#' \item{Constructor}{
#' \code{cellData(name, positive.ctr = character(0), negative.ctr = character(0),
#' expwell = character(0), norm.method =  getOption("opm.normalization.method"),
#' QC.threshold = getOption("opm.QC.threshold")}.
#' }
#' \item{Show}{signature(object="cellLoad"). Displays object content as text.}
#' \item{Accessor}{x[i]. \code{x}: a cellData object;
#' \code{i}: character, a cellData slot name.}
#' }
#' @examples
#' oneCell <- cellData(name = "Average Intensity of Nuclei",
#' positive.ctr = c("H02", "J02", "L02"),
#' negative.ctr = c("C23", "E23", "G23"))
#' oneCell
#' oneCell["name"]
#' @exportClass cellData
#' @docType class
#' @name cellData-class
#' @rdname cellData-class
#' @aliases cellData, cellData-class
setClass("cellData", slots = c(name = "character",
                         posctrwell = "character",
                         negctrwell = "character",
                            expwell = "character",
                            cellNum = "matrix",
                        origin.data = "matrix",
                          norm.data = "matrix",
                            qc.data = "data.frame",
                        norm.method = "character",
                       QC.threshold = "numeric",
                      plate.quality = "matrix",
                 plate.quality.data = "list",
                                Sig = "list"),
         prototype = list(
           QC.threshold = c(correlation = 8, zfactor = 0.5, cellnum = 50),
           norm.method = "MP"
           )
         )

#' @name cellData
#' @param name character, the analyzed item
#' @param positive.ctr a character vector, the positive control well IDs,
#' e.g. c("E05", "E06")
#' @param negative.ctr a character vector, the positive control well IDs,
#' e.g. c("B05", "B06")
#' @param expwell include all wells except control and neglect.well if NULL
#' @param norm.method character the normalization method.
#' @param QC.threshold numeric, the thresholds in the quality control.
#' @rdname cellData-class
#' @return a cellData object
#' @aliases cellData, cellData-method
#' @export
cellData <- function(name, positive.ctr = character(0),
                     negative.ctr = character(0), expwell = character(0),
                     norm.method =  getOption("opm.normalization.method"),
                     QC.threshold = getOption("opm.QC.threshold")) {
    if (any(sapply(combn(list(positive.ctr, negative.ctr, expwell),
                 2, simplify = FALSE),
           function(x) length(Reduce(intersect, x)))))
      stop("Intersection exists between different types of wells")

    name <- gsub(" ", ".", name)
    new("cellData", name = name,
        posctrwell = positive.ctr, negctrwell = negative.ctr, expwell = expwell,
        norm.method = norm.method[1], QC.threshold = QC.threshold)
}

## Getter
#' @rdname cellData-class
#' @param x a cellData object
#' @param i a requested slot name
setMethod("[", signature = c(x = "cellData", i = "character"), function(x, i) {
    if(i == "name")               { return(x@name)              } else {}
    if(i == "origin.data")        { return(x@origin.data)       } else {}
    if(i == "norm.data")          { return(x@norm.data)         } else {}
    if(i == "qc.data")            { return(x@qc.data)           } else {}
    if(i == "plate.quality")      { return(x@plate.quality)     } else {}
    if(i == "plate.quality.data") { return(x@plate.quality.data)} else {}
    if(i == "Sig")                { return(x@Sig)               } else {}
    if(i == "norm.method")        { return(x@norm.method)       } else {}
    if(i == "posctrwell")         { return(x@posctrwell)        } else {}
    if(i == "negctrwell")         { return(x@negctrwell)        } else {}
    if(i == "expwell")            { return(x@expwell)           } else {}
    if(i == "cellNum")            { return(x@cellNum)           } else {}
    if(i == "QC.threshold")       { return(x@QC.threshold)       } else {}
})
## Show method
#' @rdname cellData-class
#' @param object a cellData class
#' @aliases show,cellData-method
setMethod("show", signature = (object = "cellData"), function(object) {
    cat("An object of cellData class.\n")
    cat("Parameter: ", object@name,".\n")
    cat("Raw Data: ")
    if (length(object@origin.data)) {
        str(object@origin.data)
        cat("\n")
    } else {
        cat("NULL.\n")
    }
    cat("Data Normalization: ",
        ifelse(length(object@norm.data), "Done", "ToDo..."), ";\n", sep = "")
    cat("     method: ",
        paste(object@norm.method, collapse=" "),
        ".\n", sep = "")
    cat("Quality Control: ",
        ifelse(length(object@qc.data), "Done", "ToDo..."), ".\n", sep = "")
    cat("Significant hits: \n")
    if (length(object@Sig)) {
        show(object@Sig$stats)
    } else {
        cat("      ToDo...\n")
    }
})

#' Data importing
#' 
#' Extracts data of a specific type in a list of \code{expData} objects to
#' initialize a \code{cellData} object.
#' @param object a cellData object
#' @param lstPlates a list of expData objects
#' @param positive.ctr   a character vector, the positive control well IDs,
#' e.g. c("E05", "E06")
#' @param negative.ctr  a character vector, the positive control well IDs,
#' e.g. c("B05", "B06")
#' @param neglect.well a character vector, the neglect wells.
#' Accept regular expression, e.g. c("*02", "*23")
#' @param expwell include all wells except control and neglect.well if NULL
#' @param ... other parameters
#' @details negative.ctr accept regular expression
#' @return a \code{cellData} object, with initialized slot \code{origin.data}
#' @docType methods
#' @examples
#' data(platemap)
#' platemap$Path <- file.path(
#' system.file("Test", package = "OperaMate"), platemap$Path)
#' data(demoCell)
#' datapath <- file.path(system.file("Test", package = "OperaMate"), "Matrix")
#' lstPlates <- loadAll(cellformat = "Matrix", datapath = datapath)
#' oneCell <- cellLoad(oneCell, lstPlates, neglect.well = c("*02", "*23"))
#' str(oneCell["origin.data"])
#' @export
#' @rdname cellLoad
setGeneric("cellLoad", function(object, lstPlates, ...) {
    standardGeneric("cellLoad")
})
#' @rdname cellLoad
#' @aliases cellLoad
#' 
setMethod("cellLoad", signature = "cellData",
          function(object, lstPlates, positive.ctr = NULL, negative.ctr = NULL,
                   neglect.well = NULL, expwell = NULL) {
    
    allWell <- unique(unlist(lapply(lstPlates,
                                    function(onePlate) onePlate["wellID"])))
    
    if (!is.null(positive.ctr)) {
        object@posctrwell <- positive.ctr
    }
    if (!is.null(negative.ctr)) {
        object@negctrwell <- negative.ctr
    }
    if (is.null(expwell) & !length(object@expwell)) {
        neglect.well <-
          allWell[grep(paste(neglect.well, collapse = "|"), allWell)]
        not.expwell <- c(object@posctrwell, object@negctrwell, neglect.well)
        object@expwell <- allWell[! allWell %in% not.expwell]
    }
    allWell <- c(object@posctrwell, object@negctrwell, object@expwell)
    
    onecelldata <- matrix(NA, length(allWell), length(lstPlates))
    i <- 0
    for(onePlate in lstPlates){
        i <- i + 1
        tmpdata <- (onePlate["data"])[[object@name]]
        onecelldata[, i] <- tmpdata[match(allWell,
                                         onePlate["wellID"])]
    }
    colnames(onecelldata) <- sapply(lstPlates, function(x) x["name"])
    rownames(onecelldata) <- allWell
    object@origin.data <- onecelldata
    return(object)
})

#' Load cell number
#' 
#' @param object a cellData object
#' @param object.cellnum a cellData object for cell numbers
#' @return a \code{cellData} object, with initialized slot \code{cellNum}
#' @examples
#' data(demoCell)
#' data(demoCellNum)
#' oneCell <- cellNumLoad(oneCell, oneCellNum)
#' str(oneCell["cellNum"])
#' @docType methods
#' @export
#' @rdname cellNumLoad
setGeneric("cellNumLoad", function(object, object.cellnum) {
    standardGeneric("cellNumLoad")
})
#' @rdname cellNumLoad
#' @aliases cellNumLoad
#' 
setMethod("cellNumLoad", signature = c(object = "cellData",
                           object.cellnum = "cellData"),
          function(object, object.cellnum) {
    cell.num <- object.cellnum@origin.data
    IDRow <- match(rownames(object@origin.data), rownames(cell.num))
    IDCol <- match(colnames(object@origin.data), colnames(cell.num))
    object@cellNum <- cell.num[IDRow, IDCol]

    return(object)
})
#' Data normalization
#' 
#' Normalizes raw data based on different normalization methods.
#'
#' Method description: "MP" employes the median polish algorithm which
#' divides data by the median of their plates and wells recursively, while
#' "PMed" only divides data by the median of their plates; "Z" substracts
#' data by their plate medians, and then divides by the median absolute
#' deviations; "Ctr" divides data by the mean of their plate negative controls;
#' "None" avoids the data normalization in this step. The first three
#' methods are based on the assumption that most samples display no
#' biological effects in the assay be analyzed. They are often more
#' effective than "Ctr" method as to the high throughput screening.
#'
#' @param object a cellData object
#' @param norm.method getOption("opm.normalization.method")
#' 
#' @return a \code{celldata} object with initialized slot \code{norm.data}
#' @examples
#' data(demoCell)
#' oneCell <- cellNorm(oneCell, norm.method = "MP")
#' str(oneCell["norm.data"])
#' @export
#' @rdname cellNorm

setGeneric("cellNorm", function(object, norm.method) {
    standardGeneric("cellNorm")
})
#' @rdname cellNorm
#' @aliases cellNorm
setMethod("cellNorm", signature = "cellData",
          function(object, norm.method = getOption("opm.normalization.method")) {

    if (length(object@norm.method) > 1) {
        object@norm.method <- norm.method[1]
    } else if(length(norm.method) == 1 && object@norm.method != norm.method) {
        object@norm.method <- norm.method
    } else {}
    norm.method <- object@norm.method
    if (norm.method == "None") {
        object@norm.data <- object@origin.data
    } else if(norm.method == "PMed") {
        object@norm.data <- t( t(object@origin.data)
                              /apply(object@origin.data, 2, median, na.rm = TRUE))
    } else if(norm.method == "MP") {
        log <- capture.output(
          {tmp <- medpolish(log2(object@origin.data), na.rm = TRUE)})
        object@norm.data <- 2^(tmp$residuals)
    } else if(norm.method == "Ctr") {
        if (! is.null(object@negctrwell)) {
            lstctr <- apply(object@origin.data, 2, function(x)
                            mean(x[object@negctrwell], na.rm = TRUE))
            object@norm.data <- t(t(object@origin.data) / lstctr)
        } else {
            stop("Please specify the negative control wells")
        }
    } else if( norm.method == "Z") {
        object@norm.data <- t( t(object@origin.data)
                           - apply(object@origin.data, 2, median, na.rm = TRUE)
                              / apply(object@origin.data, 2, mad, na.rm = TRUE) )
    } else {
        stop("Set your own normalization rule to object@norm.method
             and the corresponding normalized data to object@norm.data.")
    }
    
    return(object)
})


#' Quality control
#' 
#' Checks quality of all plates and then wells.
#'
#' Requires three or more replicated samples.
#' 
#' @param object a cellData object
#' @param qcType the type of quality control
#' @param qc.threshold quality control thresholds
#' @param replace.badPlateData if TRUE,
#' replace the values of bad plate by their replicates
#' @param plot if TRUE, plot figures
#' @param outpath directory of output figures, default: getOption("opm.outpath")
#' @param ... arguments for the graphic device
#'
#' @details qcType include c("plateCorrelation", "wellSd", "zFactor", "cellNumber"),
#' An example of qc.threshold is c(correlation = 0.8, zfactor = 0.5, cellnumber = 50).
#' @return a \code{cellData} object with intialized slot \code{qc.data},
#' \code{plate.quality} and \code{plate.quality.data}.
#' @rdname cellQC
#' @examples
#' data(demoCell)
#' op <- options("device")
#' options("device" = "png")
#' oneCell <- cellQC(oneCell, qcType = c("plateCorrelation", "wellSd", "cellNumber"),
#' qc.threshold = c(correlation = 0.7), outpath = tempdir())
#' options(op)
#' str(oneCell["qc.data"])
#' str(oneCell["plate.quality"])
#' @import pheatmap
#' @export
#' 
setGeneric("cellQC",function(object, qcType = NULL, qc.threshold = NULL,
                             replace.badPlateData = TRUE,
                             plot = TRUE, outpath = getOption("opm.outpath"), ...) {
    standardGeneric("cellQC")
})
#' @rdname cellQC
#' @aliases cellQC
setMethod("cellQC", signature = "cellData",
          function(object, qcType = getOption("opm.QC.type"),
                   qc.threshold = getOption("opm.QC.threshold"),
                   replace.badPlateData = getOption("opm.replace.badPlateData"),
                   plot = TRUE, outpath = getOption("opm.outpath"), ...) {

    tmp <- getOption("opm.QC.threshold")
    qc.threshold <- c(qc.threshold, tmp[!names(tmp) %in% names(qc.threshold)])
    object@QC.threshold <- qc.threshold
          
    out <- plateQC(object, type = qcType,
                   qth = qc.threshold["correlation"],
                   zth = qc.threshold["zfactor"],
                   plot = plot, outpath = outpath, prefix = object@name, ...)
    object@plate.quality <- out$plate.quality
    object@plate.quality.data <- out$plate.quality.data
    object@qc.data <- wellQC(out$welldata, out$plate.quality, cell.num = object@cellNum,
                             type = qcType, replace.badPlateData = replace.badPlateData,
                             nth = qc.threshold["cellnumber"],
                             plot = plot, outpath = outpath, prefix = object@name, ...)
    
    return(object)
})


#' Hit identification
#' 
#' Detects samples those are most different from the negative controls.
#' 
#' @param object a cellData object
#' @param method method = c("stable","ksd","kmsd").
#' Details are referred in the vignette.
#' @param th numeric, the thresholds. It can be one threshold  for both
#' high and low expressed hit or two thresholds for each respectively.
#' @param thPval numeric, threshold of pvalues in the t-test
#' between the sample and control replicates
#' @param digits integer, the number of digits used to show the thresholds
#' @param adjust.method pvalue correction method
#' @param plot plot QQ-plot when method is "stable" if TRUE.
#' @param outpath directory of output figures, default: getOption("opm.outpath")
#' @param ... arguments of the graphic device
#' @return a \code{cellData} object with initialized slot \code{Sig}.
#' @rdname cellSig
#' @examples
#' data(demoCell)
#' op <- options("device")
#' options("device" = "png")
#' oneCell <- cellSig(oneCell, method = "stable", th = c(0.05, 0.05),
#' outpath = tempdir())
#' options(op)
#' names(oneCell["Sig"])
#' @import stats fBasics stabledist
#' @export
#' 
setGeneric("cellSig", function(object, method = c("stable", "ksd", "kmsd"),
                               th = NULL, thPval = 0.05, digits=3,
                               adjust.method = p.adjust.methods, plot = TRUE,
                               outpath = getOption("opm.outpath"), ...) {
    standardGeneric("cellSig")
})
#' @rdname cellSig
#' @aliases cellSig
setMethod("cellSig", signature = "cellData", function(
                       object, method = c("stable","ksd","kmsd"),
                       th = NULL, thPval = 0.05, digits = 3,
                       adjust.method = p.adjust.methods, plot = TRUE,
                       outpath = getOption("opm.outpath"), ...) {
    
    if (is.null(th) & method == "stable") {
        th <- c(0.05, 0.05)
    } else if (is.null(th) & method != "stable") {
        th <- c(3, 3)
    } else if (length(th) == 1) {
        th <- c(th, th)
    } else {
        th <- th[1:2]
    }
    names(th) <- c("thLow", "thHigh")

    if (is.logical(object@qc.data[,"pass.wellQC"])) {
        goodWellID <- which(object@qc.data[,"pass.wellQC"])
    } else {
        goodWellID <- 1:nrow(object@qc.data)
    }
    repdata <- object@qc.data[goodWellID,
                              - grep("pass", colnames(object@qc.data)),
                              drop = FALSE]
    if (nrow(repdata) == 0) {
        return(object)
    }
    repdata <- repdata[, -ncol(repdata), drop = FALSE]

    if (ncol(repdata) < 3) {
        message("Cannot obtain p value with replicates less than three.")
    }
    if (is.null(thPval) | ncol(repdata) < 3) {
        flag <- 0
        well.pvalue <- NULL
    } else {
        flag <- 1
        sigP <- sigByTTest(repdata, object@negctrwell,
                           pth = thPval, adjust.method = adjust.method)
        well.pvalue <- sigP$pvalue[match(rownames(object@qc.data),
                                         rownames(sigP$sig.matrix))]
        names(well.pvalue) <- rownames(object@qc.data)
    }
    
    gene.exp <- object@qc.data[goodWellID, "mean"]
    names(gene.exp) <- rownames(object@qc.data[goodWellID,])
    sig <- sigDetect(gene.exp, th["thLow"], th["thHigh"], method, 
                     digits = digits, plot = plot,
                     prefix = object@name, outpath = outpath, ...)
    threshold <- c(th, sig$threshold)
    if(flag){
        SigMat <- sig$sig.matrix & sigP$sig.matrix
    }else{
        SigMat <- sig$sig.matrix
    }
    stats <- apply(SigMat,2, sum, na.rm = TRUE)
    SigMat <- SigMat[match(rownames(object@qc.data),rownames(SigMat)), ,
                     drop = FALSE]
    rownames(SigMat) <- rownames(object@qc.data)
    
    object@Sig <- list(SigMat = SigMat,
                       threshold = threshold,
                       stats = stats,
                       pvalue = well.pvalue
                       )
    return(object)
})


#' Mean of two cellData objects
#'
#' Merges the intensities in nucleus and cytoplasm to their averages for
#' signature detection.
#'
#' @param name the name of mean cellData object
#' @param cell1  one cellData object
#' @param cell2 another celldata object
#'
#' @return the mean cellData object
#' @export
#' @rdname cellMean
#' @examples
#' data(demoCell)
#' meanCell <- cellMean(oneCell, oneCell, "meanCell")
#' meanCell
setGeneric("cellMean", function(cell1, cell2, name) {
    standardGeneric("cellMean")
})
#' @rdname cellMean
#' @aliases cellMean
setMethod("cellMean",
signature = c(cell1 = "cellData", cell2 = "cellData", name = "character"),
function(cell1, cell2, name){
    oneCell <- cellData(name = name,
                     expwell = cell1@expwell,
                positive.ctr = cell1@posctrwell,
                negative.ctr = cell1@negctrwell,
                QC.threshold =
                        apply(rbind(cell1@QC.threshold, cell2@QC.threshold),
                              1, max, na.rm = TRUE),
                 norm.method =
                        unique(c(cell1@norm.method, cell2@norm.method))
                        )

    oneCell@origin.data <- 
      (cell1@origin.data + cell2@origin.data) / 2
    oneCell@norm.data <-
      (cell1@norm.data + cell2@norm.data) / 2
    oneCell@plate.quality <-  ifelse(
      (cell1@plate.quality & cell2@plate.quality), TRUE, FALSE)
    oneCell@plate.quality.data <- c(cell1@plate.quality.data,
                                    cell2@plate.quality.data)
    oneCell@cellNum <- (cell1@cellNum + cell2@cellNum) / 2
    ID <- grep("pass", colnames(cell1@qc.data))
    if(is.logical(cell1@qc.data[, ID[length(ID)]]) &
       is.logical(cell2@qc.data[, ID[length(ID)]])){
        qc.pass <- cell1@qc.data[, ID] & cell2@qc.data[, ID]
    } else{
        qc.pass <- cell1@qc.data[, ID[-length(ID)]] &
          cell2@qc.data[, ID[-length(ID)]]
        qc.pass <- data.frame(qc.pass, pass.wellQC = rep("-", nrow(qc.pass)))
    }
    colnames(qc.pass) <- colnames(cell1@qc.data[, ID])
    oneCell@qc.data <- data.frame(
      (cell1@qc.data[, -ID] + cell2@qc.data[, -ID]) / 2,
      qc.pass)
    
    return(oneCell)
})




