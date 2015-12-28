#' @name demoData
#' @title Examples of tables and cellData objects
#' @docType data
#' @rdname data
#' @keywords data
NULL
#' @name platemap
#' @aliases platemap
#' @rdname data
#' @section platemap:
#' \describe{
#' \item{Description}{The experiment information of each Columbus
#' analysis report. This table is required only if the report formats are not
#' standarded. See \code{\link{loadAll}} for more information.}
#'\item{Format}{data.frame with the following required column names:
#' \describe{
#' \item{\code{FileName}:}{character, the name of the report.}
#' \item{\code{Format}:}{character, only ''Tab'' and ''Matrix'' are supported in
#' the current version.}
#' \item{\code{Barcode}:}{character, the barcode of the plates.}
#' \item{\code{RepID}:}{character, the ID to distinguish the replicated plates.}
#' \item{\code{Path}:}{character, the full path of the report.}
#' }
#' }
#' }
#' @examples
#' data(platemap)
#' str(platemap)
#' @return platemap: a data frame
NULL

#' @name oneCell
#' @aliases oneCell
#' @aliases demoCell
#' @rdname data
#' @section oneCell:
#' \describe{
#' \item{Description}{
#' \code{oneCell} is a \code{cellData} object used in the examples of the
#' package.}
#'}
#' @examples
#' data(demoCell)
#' oneCell
#' @return oneCell: a \code{cellData} object
NULL
#'
#' oneCellNum
#' @name onCellNum
#' @aliases oneCellNum
#' @rdname data
#' @section oneCellNum:
#' \describe{
#' \item{Description}{
#' \code{oneCellNum} is a \code{cellData} object storing the cell numbers.}
#'}
#' @examples
#' data(demoCellNum)
#' @return oneCellNum: a \code{cellData} object
NULL

## Quality control of each sample
##
## Performs quality control to each well based on its duplicated samples.
## An internal function of \code{\link{cellQC}}.
## @return a matrix containing the data after quality control, with the rows are
## "barcode:wellID" and columns are the data of each replicated samples and
## their means and if they have passed the quality control.
##
wellQC <- function(welldata, plate.quality, cell.num = NULL,
                   type = getOption("opm.QC.type"),
                   replace.badPlateData = TRUE, nth = 50,
                   plot = TRUE, outpath = getOption("opm.outpath"),
                   prefix = "", ...) {
    welldata.new <- welldata
    type <- c("wellSd", "cellNumber") %in% type
     if (type[2] & length(cell.num) == 0) {
         warning("Cell numbers are not provided.
                  Quality control will not include cell number filtering.")
         type[2] <- FALSE
     }

     barcode <- sapply(strsplit(rownames(welldata.new), ":"), function(x) x[1])
     pass.plateQC <- plate.quality[match(barcode, rownames(plate.quality)), ]
     rownames(pass.plateQC) <- rownames(welldata.new)

     if(replace.badPlateData) {
         ID <- rownames(plate.quality)[!(apply(plate.quality, 1, all))]
         for (i in ID) {
             rep.status <- plate.quality[i, ]
             if (all(!rep.status)) {
                 welldata.new[barcode == i, ] <- NA
             } else {
                 welldata.new[barcode == i, !rep.status] <-
                   apply(welldata.new[barcode == i, rep.status, drop = FALSE],
                         1, mean, na.rm = TRUE)
             }
         }
     }
     if (type[2]) {
         mat.cellnum <- cellNumberReshape(cell.num)
         pass.cellNum <- (mat.cellnum > nth)
         welldata.new[!pass.cellNum] <- NA
         if(replace.badPlateData) {
             welldata.new <- t(apply(welldata.new, 1, function(x) {
                 x[is.na(x)] <- mean(x, na.rm = TRUE)
                 x
             }))
         }
     }

     mean.data <- apply(welldata.new, 1, mean, na.rm = TRUE)
     if (type[1] & ncol(welldata) < 3) {
         warning("Only less than 3 replicates imported.
                  Quality control will not include wellSd method.")
     }
    if(all(is.nan(na.omit(mean.data)))) {
        warning("All plates fail quality control.")
        type[1] <- FALSE
    }
    if (!type[1] | (type[1] & ncol(welldata) < 3)) {
        gene.quality <- data.frame(welldata, mean = mean.data,
                                   pass = pass.plateQC,
                                   pass.wellQC = rep("-", length(mean.data)),
                                   stringsAsFactors = FALSE)
        return(gene.quality)
    }
    sd.data <- apply(welldata.new, 1, sd, na.rm = TRUE)
    useID <- (!is.na(mean.data)) & (!is.na(sd.data))
    fit.mean <- ifelse(mean.data[useID] > 1,
                       mean.data[useID], 1 / mean.data[useID]) 
    fit.sd <- sd.data[useID]

    tmp <- split(fit.sd,
                 cut(fit.mean, breaks = seq(1, max(fit.mean), length.out = 50)))
    IDs <- (sapply(tmp, length) < 10)
    while (any(IDs)) {
        id <- which(IDs)[1]
        if (id != 1) {
            tmp[[id - 1]] <- c(tmp[[id - 1]], tmp[[id]])
        } else {
            tmp[[2]] <- c(tmp[[1]], tmp[[2]])
        }
        IDs <- IDs[-id]
        tmp <- tmp[-id]
    }
    badwell <- unlist(lapply(tmp, function(x) {
        boxplot.x <- boxplot(x, plot = FALSE)
        intersect(names(which(x > mean(x, na.rm = TRUE))), names(boxplot.x$out))
    }))
    pass.wellQC <- rep(TRUE, nrow(welldata))
    names(pass.wellQC) <- rownames(welldata)
    pass.wellQC[badwell] <- FALSE
    gene.quality <- data.frame(welldata, mean = mean.data,
                               pass = pass.plateQC,
                               pass.wellQC = pass.wellQC,
                               stringsAsFactors = FALSE)
    if (plot) {
        file <- paste(prefix, "wellQualityControl", sep = ".")
        if.close.device <- gDevice.new(file, outpath, ...)
        color <- rep("gray", length(fit.mean))
        color[!pass.wellQC] <- "red"
        boxplot(tmp)
        if(if.close.device) dev.off()
    }

    return(gene.quality)
}

## Quality control of duplicated plates
##
## Performs quality control based on duplicated plates.
## An internal function of \code{\link{cellQC}}.
## @param object a cellData object
## @param type quanlity control methods c("plateCorrelation", "zFactor")
## @param qth numeric, the theshold of correlation
## @param nth numeric, the threshold of cell number
## @return a list with three components: well data matrix, quality result matrix,
## and a list of the correlation of replicates and zfactors
plateQC <- function(object, type = getOption("opm.QC.type"),
                    qth = 0.8, zth = 0.5,
                    plot = TRUE, outpath = getOption("opm.outpath"),
                    prefix = "", ...) {

    data <- object@norm.data
    pos.control <- object@posctrwell
    neg.control <- object@negctrwell

    type <- c("plateCorrelation", "zFactor") %in% type
    
    if (is.null(data)) {
        stop("Quality control should be performed after normalization.")
    }
    if (type[2] & !(length(pos.control) & length(neg.control))) {
        warning("Zfactors require positive and negative controls.")
        type[2] <- FALSE
    }
    split.rep <- split(colnames(data),
                       unlist(lapply(colnames(data),
                                     function(x) {unlist(strsplit(x, "-"))[1]}))
                       )
    lst.rep.data <- lapply(names(split.rep), function(barcode) { 
        rep.data <- data[, split.rep[[barcode]], drop = FALSE]
        colnames(rep.data) <- paste0("s", seq(1, ncol(rep.data)))
        rownames(rep.data) <- paste(barcode, rownames(rep.data), sep = ":")
        rep.data
    })
    names(lst.rep.data) <- names(split.rep)
    
    if (type[1]) {
        re.plate.cortest <- plateCorrelationTest(lst.rep.data, qth = qth)
        lst.rep.cor <- re.plate.cortest$lst.rep.cor
        lst.pass.corTest <- re.plate.cortest$lst.pass.corTest
        if (plot) {
             plateCorrelationPlot(lst.rep.cor, prefix = prefix,
                                  outpath = outpath, ...)
         }
    } else {
        lst.rep.cor <- list()
        lst.pass.corTest <- list()
    }
    if (type[2]) {
        re.plate.ztest <- zFactorTest(data, pos.control, neg.control, zth = zth)
        z.factor.test <- re.plate.ztest$pass.zTest
        lst.pass.zfactor <- lapply(names(split.rep), function(barcode) { 
            rep.zfactor <- z.factor.test[split.rep[[barcode]]]
            names(rep.zfactor) <- paste0("plateQC.s", seq(1, length(rep.zfactor)))
            rep.zfactor
        })
        z.factor <- re.plate.ztest$zfactors
    } else {
        lst.pass.zfactor <- list()
        z.factor <- NULL
    }
    check.repNo <- sapply(split.rep, length)
    repNo <- max(check.repNo)
    if (length(unique(check.repNo)) > 1) {
        warning("Exists different numbers of replicates.")
        
        lst.rep.data <- addNatoSudoRep(lst.rep.data, check.repNo, repNo)
        lst.pass.corTest <- addNatoSudoRep(lst.pass.corTest, check.repNo, repNo)
        lst.pass.zfactor <- addNatoSudoRep(lst.pass.zfactor, check.repNo, repNo)
    }

    welldata <- do.call(rbind, lst.rep.data)
    if (type[1]) {
        pass.correlation <- do.call(rbind,lst.pass.corTest)
        rownames(pass.correlation) <- names(lst.pass.corTest)
        plate.quality <- pass.correlation
    } else {
        pass.correlation <- NULL
    }
    if (type[2]) {
        pass.zfactor <- do.call(rbind,lst.pass.zfactor)
        rownames(pass.zfactor) <- names(lst.pass.zfactor)
        plate.quality <- pass.zfactor
    } else {
        pass.zfactor <- NULL
    }
    if (sum(type) == 2) {
        plate.quality <- pass.correlation & pass.zfactor
    } else if (sum(type) == 0) {
        plate.quality <- matrix(TRUE, nrow = length(lst.rep.data), ncol = repNo)
        rownames(plate.quality) <- names(lst.rep.data)
        colnames(plate.quality) <-  paste0("plateQC.s", seq(1, repNo))
    }
    plate.quality.data <- list(plate.correlation = lst.rep.cor,
                               plate.zfactor = z.factor)
    return(list(welldata = welldata, plate.quality = plate.quality,
                plate.quality.data = plate.quality.data))
}

## Hit detection methods
##
## An internal function of \code{\link{cellSig}}.
## @param gene.exp vector, the data of each sample
## @pram method "stable", "ksd" or "kmsd"
## @param thLow/thHigh numeric, the thresholds
## @param digits integer, the number of digits used to show the results
## @param outpath a character string naming the location the figure to
## generate to
## @param prefix the title of the figure
## @param ... arguments of the graphic device
## @details ksd and kmsd select hit with values beyond the range
## [mean - kLow * sd, mean + kHigh * sd] and
## [median - kLow * mad, median + kHigh * mad], respectively. 
## @return a list including the hit matrix and the thresholds
## 
sigDetect <- function(gene.exp, thLow, thHigh,
                      method = c("stable", "ksd", "kmsd"), digits = 3,
                      plot = TRUE, prefix = "",
                      outpath = getOption("opm.outpath"), ...) {
    if (method == "ksd") {
        m.exp <- mean(gene.exp, na.rm = TRUE)
        sd.exp <-  sd(gene.exp, na.rm = TRUE)
        ths <- c(thLow = m.exp - thLow * sd.exp,
                 thHigh = m.exp + thHigh * sd.exp)
    } else if (method == "kmsd") {
        md.exp <- median(gene.exp, na.rm = TRUE)
        msd.exp <-   mad(gene.exp, na.rm = TRUE)
        ths  <- c (thLow = md.exp - thLow * msd.exp,
                   thHigh = md.exp + thHigh * msd.exp)
    } else if (method == "stable") {
        r <- fitStable(gene.exp, plot, prefix, outpath, ...) 
        ths <- c(thLow = qstable(thLow, r$estimate[1], r$estimate[2],
                   r$estimate[3], r$estimate[4]),
                 thHigh = qstable(1-thHigh, r$estimate[1], r$estimate[2],
                   r$estimate[3], r$estimate[4]))
    }

    sig.matrix <- cbind(gene.exp < ths[1], gene.exp > ths[2])
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- names(gene.exp)
    colnames(sig.matrix) <- c("Low","High")
    ths <- round(ths, digits = digits)
    names(ths) <- c("Low","High")
    return(list(sig.matrix = sig.matrix, threshold = ths))
}

## Significant samples detection by multiple t test
##
## An internal function of \code{\link{cellSig}}.
##
## Multiple t-test between samples and controls is the compulsory step for
## \code{cellSig}. It is combined with other method for hit identification.
## @param repdata matrix, the data of replicated samples
## @param ctrwell vector of character specifying the control well IDs
## @param pth numeric, the threshold of pvalue of t-test
## @param adjust.method pvalue correction method.
## @return including the hit matrix and the pvalues of all samples
## 
sigByTTest <- function(repdata, ctrwell,
                  pth = 0.05, adjust.method = p.adjust.methods) {

    g.ctrexp <- na.omit(as.vector(unlist(
      repdata[grep(paste(ctrwell, collapse = "|"), rownames(repdata)), ,
              drop = FALSE])))
    lst.plate.info <- split(repdata, 
                            sapply(strsplit(rownames(repdata), ":"),
                                   function(x) x[1])
                            )

    pval.well <- unlist(lapply(lst.plate.info, function(pmat) {
        ctrexp <-
          pmat[grep(paste(ctrwell, collapse = "|"), rownames(pmat)), ,
               drop = FALSE]
        if (nrow(ctrexp) == 0 ) {
            ctrexp <- sample(g.ctrexp, ncol(repdata)*5)
        } else  {
            ctrexp <- as.vector(unlist(ctrexp))
            ctrexp[is.na(ctrexp)] <- sample(g.ctrexp, sum(is.na(ctrexp)))
        }
        pval.well <- apply(pmat, 1, function(exp) {
            ifelse(any(is.na(exp)), NA,
                   t.test(exp, ctrexp, alternative = "less")$p.value) })
        names(pval.well) <- sapply(strsplit(names(pval.well), "[:]"),
                                   function(x) x[2])
        pval.well
    }))
    names(pval.well) <- sub("[.]", ":", names(pval.well))
    pval.well.twoside <- ifelse(pval.well < 0.5, pval.well, 1 - pval.well) * 2
    padjust.well <- p.adjust(pval.well.twoside, adjust.method)
    
    sig.ID <- (padjust.well < pth)
    sig.matrix <- cbind(pval.well < 0.5 & sig.ID,
                        pval.well > 0.5 & sig.ID)
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- rownames(repdata)
    colnames(sig.matrix) <- c("Low", "High")
    
    return(list(sig.matrix = sig.matrix, pvalue = padjust.well))
}
## Fit data to stable distribution
## 
## Fits data to a student's t distribution.
## @param gene.exp vector
## @param plot logic, if to plot the figure
## @param prefix character, the name of the figure
## @param outpath a character string naming the location the figure to
## generate to. 
## @param ... arguments of the graphics device
## @param xlab character, X axis label
## @param ylab character, Y axis label
fitStable <- function(gene.exp, plot = TRUE, prefix = "",
                      outpath = getOption("opm.outpath"),
                      xlab = "Data", ylab = "Stable Distribution", ...) {
    gene.exp <- na.omit(gene.exp)
    r <- stableFit(gene.exp, gamma = log(max(gene.exp)),
                   delta = median(gene.exp), doplot = FALSE)@fit
    if (plot) {
        file <- paste(prefix, "dataFitting", sep = ".")
        if.close.device <- gDevice.new(file, outpath, ...)
        qqplot(qstable(ppoints(100), r$estimate[1],
                       r$estimate[2], r$estimate[3],r$estimate[4]),
               sort(gene.exp),
               xlab = xlab, ylab = ylab,
               main = "QQ plot with stable distribution"
               )
        abline(0,1)
        if (if.close.device) dev.off()
    }
    return(r)
}

## Visualize data of one plate
##
## An internal function of \code{\link{plateViz}}.
## 
## @param data matrix, data of a plate
## @param pos.control, chracter, positive control
## @param neg.control, chracter, negative control
## @param outpath, the location of figure
## @param prefix, file name
## @param ... other auguments
## @return Invisibly a list of component by \code{pheatmap}.
## 
## @import pheatmap
## 
showPlate <- function(data, pos.control = NULL, neg.control = NULL,
                      outpath = getOption("opm.outpath"),
                      prefix = "", ...){
    file <- paste(prefix, "Visualization", sep = ".")
    if.close.device <- gDevice.new(file, outpath, ...)
    g <- pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, ...)
    if (if.close.device) dev.off()
    invisible(g)
}

#' Report generation
#'
#' Summarizes all results in the list of \code{cellData} objects,
#' and writes out a report to file.
#'
#' This function summarizes the information from all \code{cellData} objects,
#' and visualizes the number of the hists if required.
#' @param lstCells a list of cellData objects
#' @param genemap a data frame, the well-gene specification table
#' @param verbose logical, detailed data will be provided if TRUE
#' @param file the path of the file to generate to
#' @param outpath a character string naming the location the figures to
#' generate to
#' @param plot if TRUE, plot barplot
#' @param ... arguments of the graphic device
#' @return a data frame with annotated information of each well
#' @import ggplot2
#' @export
#' @examples
#' data(demoCell)
#' genemap <- read.csv(file.path(system.file("demoData", package = "OperaMate"),
#' "genemap.csv"), stringsAsFactors = FALSE)
#' report <- generateReport(list(oneCell), genemap, verbose = FALSE,
#' plot = FALSE)
#' str(report)
#' 
generateReport <- function(lstCells, genemap= NULL, verbose = FALSE, file = NULL,
                           outpath = getOption("opm.outpath"), plot = TRUE, ...) {
    if (!is.null(genemap)) {
        if(any( !c("Barcode","Well") %in% colnames(genemap))){
            warning("The mapping table should contain Barcode and Well of each gene.")
            genemap <- NULL
        }
        genemap <- genemap[genemap$Well!="",]
        matchID <- paste(genemap$Barcode, genemap$Well, sep = ":")
    }
    if (is.null(genemap)) {
        matchID <- unique(unlist(sapply(lstCells, function(cell) {
            rownames(cell["qc.data"])})))
    }
    if (!verbose) {
        meanVal <- do.call(cbind, lapply(lstCells, function(cell) {
            mean.each <-
              cell["qc.data"][, grep("mean", colnames(cell["qc.data"])), drop = FALSE]
            mean.each[match(matchID, rownames(mean.each)), ]
        }))
        pass.quality <- do.call(cbind, lapply(lstCells, function(cell) {
            q.each <- ifelse(cell["qc.data"][, "pass.wellQC", drop = TRUE],
                             "Pass", "Fail")
            names(q.each) <- rownames(cell["qc.data"])
            q.each[match(matchID, names(q.each))]
        }))
        hits <- do.call(cbind, lapply(lstCells, function(cell) {
            sigmat <- (cell["Sig"])$SigMat
            sigvec <- rep("notHit", nrow(sigmat))
            sigvec[sigmat[, 1]] <- "Low"
            sigvec[sigmat[, 2]] <- "High"
            names(sigvec) <- rownames(sigmat)
            sigvec[match(matchID, names(sigvec))]
        }))
        report <- data.frame(genemap, mean = meanVal,
                             pass.quality = pass.quality, Hits = hits)
    } else {
        cellNumID <- which(sapply(lstCells, function(x) ! is.null(x@cellNum)))
        cellnum <- cellNumberReshape(lstCells[[cellNumID[1]]]["cellNum"])
        cellnum <- cellnum[match(matchID, rownames(cellnum)), ]

        lstqc <- lapply(lstCells, function(cell) {
            qc.data <- cell["qc.data"]
            colnames(qc.data) <- paste(cell["name"], colnames(qc.data), sep = ".")
            qc.data <- qc.data[match(matchID, rownames(qc.data)), , drop = FALSE]
        })
        qc.data <- Reduce(cbind, lstqc)
        
        lstsigdata <- lapply(lstCells, function(cell) {
            sig.data <- cell["Sig"]$SigMat
            colnames(sig.data) <- paste(cell["name"], colnames(sig.data), sep = ".")
            sig.data[match(matchID, rownames(sig.data)), , drop = FALSE]
        })
        sig.data <- Reduce(cbind, lstsigdata)

        report <- data.frame(genemap, qc.data, sig.data)
    }
    threshold <- unlist(lapply(lstCells, function(obj) {
        th <- obj["Sig"]$threshold[c("Low", "High")]
        return(th)
    }))
    stats <- unlist(lapply(lstCells, function(obj) {
        st <- obj["Sig"]$stats[c("Low","High")]
        return(st)
    }))
    statMat <- cbind(c("threshold", "number"), rbind(threshold, stats))
    if (!is.null(file)) {
        write("Analysis Report", file = file.path(outpath, file))
        suppressWarnings(write.table(statMat, row.names = FALSE, sep = "\t",
                                     file = file.path(outpath, file),
                                     quote = FALSE, append = TRUE))
        cat("\n",file = file.path(outpath, file), append = TRUE)
        suppressWarnings(write.table(report, file = file.path(outpath, file),
                                     append = TRUE,  sep = "\t", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }
    if (plot) {
        if.close.device <- gDevice.new("Statistics", outpath, ...)
        count <- NULL; subType <- NULL; Type <- NULL
        ## Just for avoiding warnings in Rcheck
        stats.ggplot <- data.frame(count = stats,
                                   subType = names(stats),
                                   Type = rep(sapply(lstCells,
                                     function(obj) obj["name"]), each = 2)
                                   )
        xlab.ggplot <- apply(stats.ggplot, 1, function(x)
                             sub(paste0(x["Type"], "."), "", x["subType"]))
        print(ggplot(data = stats.ggplot, aes(x = subType, y = count)) +
              geom_bar(aes(fill = Type), stat = "identity") +
              scale_x_discrete(labels = xlab.ggplot ) +
              geom_text(aes(label = count)) +
              xlab(""))
        if (if.close.device) dev.off()
    }
    return(report)
}

## Inner function
gDevice.new <- function(filename, outpath = getOption("opm.outpath"), ...) {
    gDevice <- getOption("device")
    if.close.device <- FALSE
    if(is.character(gDevice)) {
        gDevice <- if(exists(gDevice, .GlobalEnv)) get(gDevice, .GlobalEnv)
        else if(exists(gDevice, asNamespace("grDevices")))
          get(gDevice, asNamespace("grDevices"))
        else stop(gettextf("device '%s' not found", gDevice), domain=NA)
    }
    pars <- list(...)
    device.pars <- formals(gDevice)
    pars <- pars[names(pars) %in% names(device.pars)]
    filename <- file.path(outpath,
                          filename)
    if(any(names(device.pars) == "title")) {
        do.call(dev.new, c(title = filename, pars))
    } else if(any(names(formals(gDevice)) == "filename")) {
        if.close.device <- TRUE
        ext <- tools::file_ext(formals(gDevice)$filename)
        do.call(dev.new, c(filename = paste(filename, ext, sep="."), pars))
    } else {
        do.call(dev.new, pars)
    }
    return(if.close.device)
}

#' Plate information extraction
#' 
#' Extract plate information from file names.
#'
#' @param vec.files a vector of file names
#' @param egFilename a file name example
#' @details An example of egFilename = list(eg.filename = "0205-s2-01.txt",
#' rep.id = "s2", exp.id = "01", sep = "-", barcode = "DSIMGA01").
#' @return a data frame of PlateID, RepID, and Barcode
#' @export
nameParser <- function(vec.files, egFilename){
    removeExt <- function (x) {
   # Adopted from limma package
        x <- as.character(x)
        n <- length(x)
        if (length(grep("\\.", x)) < n)
          return(x)
        ext <- sub("(.*)\\.(.*)$", "\\2", x)
        if (all(ext[1] == ext))
          return(sub("(.*)\\.(.*)$", "\\1", x))
        else return(x)
    }
    
    stopifnot(length(vec.files) != 0)
    eg.filename <- egFilename[["eg.filename"]]
    rep.id <- egFilename[["rep.id"]]
    exp.id <- egFilename[["exp.id"]] 
    sep <- egFilename[["sep"]]
    barcode <- egFilename[["barcode"]]
    
    eg.file.elements <- unlist(strsplit(removeExt(eg.filename), sep))
    plateID.position <- match(exp.id, eg.file.elements)
    repID.position <- match(rep.id, eg.file.elements)

    filename.content <- removeExt(basename(vec.files))
    tmp <- strsplit(filename.content, sep)
    vec.plate.id <- sapply(tmp, function(x) x[plateID.position])
    vec.rep.id <- sapply(tmp, function(x) x[repID.position])

    ## Deal with Barcode01 vs filename s1-1.txt problem
    if (length(grep("^[0-9]+$", exp.id)) == 0) { ## expID character
        exp.id <- as.character(exp.id)
        newBarcode <-
          sapply(vec.plate.id, function(x) {sub(exp.id, x, barcode) })
    } else {
        m <- gregexpr('0*[0-9]+$', barcode)
        egFPlateID <- regmatches(barcode, m)[[1]]
        if (length(egFPlateID) & length(grep(exp.id, egFPlateID))){
            if(egFPlateID == exp.id){
                newBarcode <- sapply(vec.plate.id, function(x){
                    sub(egFPlateID, x, barcode) })
            }else{
                vec.plate.id <- as.numeric(vec.plate.id)
                FormatPlateID <- formatC(vec.plate.id, width = nchar(egFPlateID),
                                         flag = "0")
                newBarcode <- sapply(FormatPlateID, function(x){
                    sub(egFPlateID, x, barcode) })
            }
        } else {
            stop("The expID and Barcode does not match!")
        }
    }

    return(data.frame(PlateID = vec.plate.id, RepID = vec.rep.id,
                      Barcode = newBarcode, row.names = NULL, check.names = FALSE))
}

## Calculate z factors of each screen
zFactorTest <- function(data, pos.control, neg.control, zth = 0.5) {
    if (any(! c(pos.control, neg.control) %in% rownames(data))) {
        stop("Control wells are not provided in the data.")
    }
    if (length(intersect(pos.control, neg.control)) > 0) {
        stop("Wrong control wells input.")
    }
    z.factor <- apply(data, 2, function(x) {
        mup <- mean(x[pos.control], na.rm = TRUE)
        mun <- mean(x[neg.control], na.rm = TRUE)
        sigmap <- sd(x[pos.control], na.rm = TRUE)
        sigman <- sd(x[neg.control], na.rm = TRUE)
        z <- 1 - 3 * (sigmap + sigman) / abs(mup - mun)
        z
    })
    pass.zTest <- (z.factor > zth)
    z.factor.info <- list(pass.zTest = pass.zTest, zfactors = z.factor)
    return(z.factor.info)
}

cellNumberReshape <- function(cellnum) {
    split.rep <- split(colnames(cellnum),
                       unlist(lapply(colnames(cellnum),
                                     function(x) {unlist(strsplit(x, "-"))[1]}))
                       )
    lst.rep.cellnum <- lapply(names(split.rep), function(barcode) {
        rep.data <- cellnum[, split.rep[[barcode]], drop = FALSE]
        colnames(rep.data) <- paste0("s", seq(1, ncol(rep.data)))
        rownames(rep.data) <- paste(barcode, rownames(rep.data), sep = ":")
        rep.data
    })
    names(lst.rep.cellnum) <- names(split.rep)
    check.repNo <- sapply(split.rep, length)
    if (length(unique(check.repNo)) > 1) {
        warning("Exists different numbers of replicates.")
        repNo <- max(check.repNo)
        lst.pass.cellNum <- addNatoSudoRep(lst.pass.cellNum, check.repNo, repNo)
    }
    cellNum <- do.call(rbind, lst.rep.cellnum)
    return(cellNum)
}

plateCorrelationTest <- function(lst.rep.data, qth = 0.8) {
    lst.rep.cor <- lapply(lst.rep.data, function(rep.data) {
        if (ncol(rep.data) == 1) {
            rep.cor <- as.matrix(1)
        } else {
            rep.cor <- cor(rep.data, use = "na.or.complete", method = "pearson")
            }
        rep.cor
    })
    lst.pass.corTest <- lapply(lst.rep.cor, function(rep.cor) {
        pass.plateQC <- c(plateQC.s = rep(TRUE, nrow(rep.cor)))
        rep.cor.pass <- (rep.cor > qth)
        while (any(! rep.cor.pass)) {
            pass.ID <- apply(rep.cor.pass, 1, all)
            delID <-  which.min(apply(rep.cor.pass[!pass.ID, ], 1, mean))
            pass.plateQC[paste0("plateQC.",names(delID))] <- FALSE
            rep.cor.pass <- rep.cor.pass[-delID, -delID, drop = FALSE]
            if (length(rep.cor.pass) == 1) {
                pass.plateQC[paste0("plateQC.",rownames(rep.cor.pass)[1])] <- FALSE
                break
            }
        }
        return(pass.plateQC)
        })
    return(list(lst.rep.cor = lst.rep.cor, lst.pass.corTest = lst.pass.corTest))
}

plateCorrelationPlot <- function(lst.rep.cor, prefix = "",
                                 outpath = getOption("opm.outpath"), ...){
    if (length(lst.rep.cor) ==0 )
      return()
    check.repNo <- sapply(lst.rep.cor, nrow)
    repNo <- unique(check.repNo)
    if (length(repNo) > 1) {
        warning("Exists different numbers of replicates.")
        repNo <- max(repNo)
        lst.rep.cor[check.repNo != repNo] <-
          lapply(lst.rep.cor[check.repNo != repNo],
                 function(rep.cor){
                     rep.cor.new <- matrix(NA, repNo, repNo)
                     rep.cor.new[1:nrow(rep.cor), 1:ncol(rep.cor)] <- rep.cor
                     rep.cor.new
                 })
    }

    cor.data <- do.call(cbind, lst.rep.cor)
    col.labels <- rep("", ncol(cor.data))
    col.labels[seq(ceiling(repNo/2), ncol(cor.data), repNo)] <-
      names(lst.rep.cor)
    
    file <- paste(prefix, "plateQualityControl", sep = ".")
    if.close.device <- gDevice.new(file, outpath, ...)
    g <- pheatmap(cor.data, color = getOption("opm.heatmap.color"),
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  legend_labels = "Pearson correlation",
                  show_rownames = TRUE, show_colnames = TRUE,
                  main = "Pearson correlation between replicates",
                  fontsize = getOption("opm.heatmap.fontsize"),
                  fontsize_row = getOption("opm.heatmap.fontsize_row"),
                  fontsize_col = getOption("opm.heatmap.fontsize_col"),
                  gaps_col = seq(0, ncol(cor.data), repNo),
                  labels_col = col.labels
                  )
    if(if.close.device) dev.off()
    return(invisible(g))
}

addNatoSudoRep <- function(lst.data, check.repNo, max.repNo){
    if (length(lst.data) == 0) {
        return(list())
    }
    lst.data[check.repNo != max.repNo] <-
      lapply(lst.data[check.repNo != max.repNo],
             function(rep.data) {
                 if (is.null(nrow(rep.data))) {
                     rep.data.na <-
                       c(rep.data, rep(NA, max.repNo - length(rep.data)))
                 } else {
                     rep.data.na <- 
                       data.frame(rep.data, matrix(NA, nrow(rep.data),
                                                   max.repNo - ncol(rep.data)))
                 }
             rep.data.na
             })
    lst.data
}

#' The barplot of enrichment functions
#' @param chart data frame, the functional annotation chart
#' @param prefix character, the prefix of figure name
#' @param type selected domains from chart, e.g. BP.
#' @param fill color of the bars
#' @param outpath directory of output figures, default: getOption("opm.outpath")
#' @param ... other arguments for graphical devices
#' @return Invisibly the ggplot2 function for barplot
#' @docType methods
#' @examples
#' data(demoCell)
#' genemap <- read.csv(file.path(system.file("demoData", package = "OperaMate"),
#' "genemap.csv"), stringsAsFactors = FALSE)
#' chart <- cellSigAnalysis(oneCell, genemap, organism = "mmusculus")
#' op <- options("device")
#' options("device" = "png")
#' cellSigAnalysisPlot(chart, type = "BP", outpath = tempdir())
#' options(op)
#' @import ggplot2
#' @export 
cellSigAnalysisPlot <- function(chart, prefix = "", type = NULL, fill = "steelblue",
                                outpath = getOption("opm.outpath"), ...) {
    if (!is.null(type)) {
        chart <- chart[chart[, "domain"] %in% type, ]
    }
    if.close.device <- gDevice.new(paste(prefix, "HitFunctions", sep = "."),
                                   outpath, ...)
    data <- data.frame(term = chart[, "term.name"],
                     pvalue = (-log10(chart[, "p.value"])))
    term <- NULL; pvalue <- NULL
    g <- ggplot(data, aes(x = term, y = pvalue)) +
      geom_bar(stat = "identity", fill = fill) +
        theme_minimal() + coord_flip()
    print(g)
    if (if.close.device)   dev.off()

    invisible(g)
}



