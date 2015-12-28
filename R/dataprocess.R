#' Data process and analysis pipeline
#' 
#' A systematical pipeline for opera data importing, normalization, quality
#' control, hit detection, analysis, and visualization.
#' @param configFile the location of the file specifying all parameters
#' @param gDevice the graphics device 
#' @param ... addition arguments for graphics devices
#' @return a list of three components:  a list of cellData objects, the
#' annotated table of each well, and the enrichment analysis table
#' @export
#' @examples
#' configFile <- file.path(system.file("demoData", package = "OperaMate"), "demoParam.txt")
#' operaReport <- operaMate(configFile, gDevice = "png")
#' head(operaReport$report)
#' 
operaMate <- function(configFile, gDevice = "png", ...) {

    ## To avoid warning in RCheck
    platemap <- NULL; eg.filename <- NULL; exp.id <- NULL
    rep.id <- NULL;  cellnames <- NULL; barcode <- NULL; well.digits <- NULL
    genemap <- NULL; cellformat <- NULL; cell.digits <- NULL
    correlationTh <- NULL; zfactorTh <- NULL; cellnumberTh <- NULL
    positive.control <- NULL; negative.control <- NULL
    neglect.well <- NULL; case.well <- NULL; norm.method <- NULL
    opm.qcType <- NULL; opm.QC.threshold <- NULL; sep <- NULL
    if.replace.badPlateData <- NULL; sig.method <- NULL
    sig.threshold <- NULL; sig.pvalue.threshold <- NULL
    verbose <- NULL; summaryFile <- NULL; organism <- NULL
    functionFile <- NULL; egFilename <- NULL; datapath <- NULL

    configParser <- function(configFile) {
        config <- file(configFile, open = "r")
        lines <- readLines(config)
        ID <- grep("^#", lines)
        lines <- gsub(" *", "", lines[-ID])
        lstlines <- strsplit(lines, ":")
        i <- 0
        for (line in lstlines) {
            i <- i + 1
            if (length(line) == 1) {
                eval(parse(text = paste(line[1]," = NULL", sep = "")))
            } else {
                tmp <- paste("\"", unlist(strsplit(line[2], ",")), "\"", sep = "")
                if (length(tmp) > 1){
                    tmp <- paste("c(", paste(tmp, collapse = ","), ")", sep = "")
                }
                eval(parse(text = paste("assign(\"", line[1], "\", ", tmp,
                             ", parent.env(environment()))", sep = "")))
            }
        }
        if (is.null(datapath)){
            stop("You must provide the location of files!")
        } else if (datapath == "operaMateDemoLocation") {
            assign("datapath",
                   file.path(system.file("Test", package = "OperaMate"), "Matrix"),
                   parent.env(environment()))
        }
        if (!file.exists(datapath)) {
            stop("The location does not exists!")
        }
        if (is.null(outpath)) {
            assign("outpath", getOption("opm.outpath"), parent.env(environment()))
        }
        if (outpath == "operaMateDemoOutput") {
            assign("outpath", tempdir(), parent.env(environment()))
        }
        if (!is.null(platemap)) {
            assign("platemap", read.csv(platemap, stringsAsFactors = FALSE),
                   parent.env(environment()))
        }
        if (is.null(genemap)) {
            warning("No well-gene specification file provide.
                 Functional analysis will not perform.")
        } else {
            if (genemap == "operaMateDemoGenemap")
              genemap <- file.path(system.file("demoData", package = "OperaMate"),
                                   "genemap.csv")
            assign("genemap", read.csv(genemap, stringsAsFactors = FALSE),
                   parent.env(environment()))
        }

        if (is.null(eg.filename) | is.null(rep.id) | is.null(exp.id)
            | is.null(sep) | is.null(barcode)) {
            warning("File format  is not clear. The standard format will be used.")
            egFilename <- getOption("opm.filename.example")
        } else{
            egFilename <- list(eg.filename = eg.filename, rep.id = rep.id,
                               exp.id = exp.id, sep = sep, barcode = barcode)
        }
        assign("egFilename", egFilename, parent.env(environment())) 
        if (is.null(cellnames)) {
            stop ("Please specify the terms to analyze!")
        }
        opm.QC.threshold <- vector()
        if (!is.null(correlationTh))
          opm.QC.threshold["correlation"] <- correlationTh
        if (!is.null(zfactorTh))
          opm.QC.threshold["zfactor"] <- zfactorTh
        if (!is.null(cellnumberTh))
          opm.QC.threshold["cellnumber"] <- cellnumberTh
        mode(opm.QC.threshold) <- "numeric"
        assign("opm.QC.threshold", opm.QC.threshold, parent.env(environment()))
        mode(sig.threshold) <- "numeric"
        mode(sig.pvalue.threshold) <- "numeric"
        mode(well.digits) <- "numeric"
        mode(if.replace.badPlateData) <- "logical"
        mode(verbose) <- "logical"
        assign("sig.threshold", sig.threshold, parent.env(environment()))
        assign("sig.pvalue.threshold", sig.pvalue.threshold, parent.env(environment()))
        assign("well.digits", well.digits, parent.env(environment()))
        assign("if.replace.badPlateData", if.replace.badPlateData, parent.env(environment()))
        assign("verbose", verbose, parent.env(environment()))
        
        close(config)
    }

    configParser(configFile) #Parse all parameters from text file

    if (is.null(outpath)) {
        outpath <- getOption("opm.outpath")
    }
    if (!file.exists(outpath)) {
        dir.create(outpath)  #Create the output folder
    }

    op <- options("device")
    options("device" = gDevice)
    
    message("[",format(Sys.time(), "%m-%d-%Y %T"),"]")
    message(" OperaMate Data Processing & Analysis")  #Log head
    message("********************************************************")

    message("Loading data ...")
    lstPlates <- loadAll(cellformat = cellformat, datapath = datapath,
                         egFilename = egFilename, well.digits = well.digits,
                         platemap = platemap)
    
    cellnames.1 <- cellnames[cellnames %in% names(lstPlates[[1]]["data"])]
    if("Average.Total.Intensity" %in% cellnames){
        cellnames.1 <- unique(c("Average.Intensity.of.Nuclei",
                                "Average.Intensity.of.Cytoplasm", cellnames.1))
    }
    lstCells <- list()
    for(cellname in c("Cells.Analyzed", cellnames.1)){
        oneCell <- cellData(cellname)
        lstCells[[cellname]] <- cellLoad(oneCell, lstPlates,
                                         positive.ctr = positive.control,
                                         negative.ctr = negative.control,
                                         neglect.well = neglect.well,
                                              expwell = case.well)
    }
    cell.cellNum <- lstCells[["Cells.Analyzed"]]
    lstCells <- lstCells[-1]
    lstCells <- lapply(lstCells, function(cell) {
        cellNumLoad(cell, cell.cellNum)
    })
    
    message("Data normalization ...")
    lstCells <- lapply(lstCells,function(cell){
        cellNorm(cell, norm.method = norm.method)
    })
    message("Data visualization ..." )
    for(cell in lstCells){
        cellViz(cell, outpath = outpath, ...)
    }

    message("Quality control ...")
    lstCells <- lapply(lstCells,function(cell){
        cellQC(cell, qcType = opm.qcType, qc.threshold = opm.QC.threshold,
               replace.badPlateData = if.replace.badPlateData,
               outpath = outpath, ...)
    })

    if("Average.Total.Intensity" %in% cellnames){
        lstCells[["Average.Total.Intensity"]] <-
          cellMean(lstCells[["Average.Intensity.of.Nuclei"]],
                   lstCells[[ "Average.Intensity.of.Cytoplasm"]],
                   "Average.Total.Intensity")
    }
    lstCells <- lstCells[cellnames]

    message("Hit detection ...")
    lstCells <- lapply(lstCells, function(cell){
        cell <- cellSig(cell, method = sig.method,
                        th = sig.threshold, thPVal = sig.pvalue.threshold,
                        adjust.method = getOption("opm.adjust.methods"),
                        digits = getOption("opm.threshold.digits"),
                        plot = TRUE, outpath = outpath, ...)
        cellSigPlot(cell, outpath = outpath, ...)
        cell
    })

    message("Annotation ... ")
    report <- generateReport(lstCells, genemap, verbose = verbose, 
                             file = summaryFile,
                             outpath = outpath, plot = TRUE, ...)

    message("Hit analysis ...")
    if (is.null(genemap)) {
        funReport <- NULL
    } else{
        funReport <- lapply(lstCells,function(cell){
            chart <- cellSigAnalysis(cell, genemap, organism,
                            file = file.path(outpath, functionFile), ...)
            cellSigAnalysisPlot(chart, prefix = cell@name, outpath = outpath,
                                ...)
            chart
        })
    }
    return(list(lstCells = lstCells, report = report, funReport = funReport))
}
