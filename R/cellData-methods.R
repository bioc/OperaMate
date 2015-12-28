#' @import stats 
#' @import methods
NULL

#' Data visualization
#' 
#' Visualize data by heatmap or boxplot.
#' @param object a cellData object
#' @param data.type c("raw", "norm), visualizing both types by default
#' @param plot c("heatmap","boxplot")
#' @param outpath directory of output figures, default: getOption("opm.outpath")
#' @param multiplot logical, the output images are placed in one figure or not
#' @param plateID numeric or character
#' @param tag character, unique tag for one figure
#' @param ctr.excluded logical, if controls are included in the visualization
#' @param ... other arguments for graphical devices and \code{pheatmap}
#' 
#' @examples
#' 
#' data(demoCell)
#' op <- options("device")
#' options("device" = "png")
#' cellViz(oneCell, data.type = c("raw", "norm"), plateID = 1:6, outpath = tempdir())
#' cellViz(oneCell, data.type = c("raw", "norm"), plateID = 1, outpath = tempdir())
#' options(op)
#' 
#' @return Invisibly a list of the values returned by \code{pheatmap} and
#' ggplot2 function for boxplot
#'
#' @details
#' By visualizing the raw data, users can observe the batch effects as
#' a large region of distinguishing color in heatmap or
#' biased distribution by boxplots. Users can also visualize thr normalized data
#' for comparison.
#' 
#' @docType methods
#' @import ggplot2 pheatmap
#' @importFrom grid pushViewport viewport grid.newpage grid.text gpar
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' @export
#' 
cellViz <- function(object, data.type = c("raw", "norm"),
                    plot = c("heatmap","boxplot"),
                    outpath = getOption("opm.outpath"),
                    multiplot = FALSE, plateID = NULL, tag = NULL,
                    ctr.excluded = TRUE, ...) {

    if(!length(object["norm.data"]) & "norm" %in% data.type) {
        msg <- "Cannot visualize normalized data before normalization! "
        if(length(data.type) == 2) {
            warning(msg,"Visualize raw data only.")
            data.type <- "raw"
        } else {
            stop(msg)
        }
    }  # parameter checking
    
    showRow <- object["expwell"]
    showCol <- colnames(object["origin.data"])
    if (!is.null(plateID)) {
        if (is.numeric(plateID)) {
            showCol <- intersect(1:length(showCol), plateID)
        } else {
            showCol <- intersect(showCol, plateID)
        }
    }
    if (length(showCol) == 1) {
        if ("raw" %in% data.type) {
            vec.plate.data <- object["origin.data"][, showCol, drop = TRUE]
        }
        if ("norm" %in% data.type) {
            vec.plate.data <- object["norm.data"][, showCol, drop = TRUE]
        }
        tmp <-  split(vec.plate.data, substr(names(vec.plate.data), 1, 1))
        mcolname <- unique(substr(names(vec.plate.data), 2, 3))
        plate.data <- do.call(rbind, lapply(tmp,function(x) {
            return(x[match(mcolname,substr(names(x),2,3))])
        }))
        colnames(plate.data) <- mcolname
        rownames(plate.data) <- names(tmp)
        g <- showPlate(plate.data, pos.control = object["posctrwell"],
                       neg.control = object["negctrwell"], outpath = outpath,
                       prefix = paste(object["name"], showCol, sep = "-"), ...)
        return(invisible(g))
    }
    
    g1 <- NULL; g2 <- NULL    
    if("heatmap" %in% plot) {
        if (! ctr.excluded) {
            showRow <- rownames(object["origin.data"])
        }
        
        if("raw" %in% data.type){
            if.close.device <-
              gDevice.new(paste0("heatmap.raw.",object["name"], tag), outpath, ...)
            g1 <- cellHeatmap(object["origin.data"][showRow, showCol], ...)
            if(if.close.device) dev.off()
        }
        if("norm" %in% data.type){
            if.close.device <-
              gDevice.new(paste0("heatmap.norm.",object["name"], tag), outpath, ...)
            g2 <- cellHeatmap(object["norm.data"][showRow, showCol], ...)
            if(if.close.device)  dev.off()
        }
    }
    
    p1 <- NULL; p2 <- NULL; p3 <- NULL; p4 <- NULL
    if("boxplot" %in% plot) {
        if (ctr.excluded) {
            positive.ctr <- NULL
            negative.ctr <- NULL
        } else {
            postive.ctr <- object["posctrwell"]
            negative.ctr <- object["negctrwell"]
        }
        if("raw" %in% data.type) {
            p1 <- cellBoxplot(object["origin.data"][showRow, showCol], "Plate",
                              positive.ctr, negative.ctr, ...)
            p2 <- cellBoxplot(object["origin.data"][showRow, showCol], "Well",
                              positive.ctr, negative.ctr, ...)
        }
        if("norm" %in% data.type) {
            p3 <- cellBoxplot(object["norm.data"][showRow, showCol], "Plate",
                              positive.ctr, negative.ctr, ...)
            p4 <- cellBoxplot(object["norm.data"][showRow, showCol], "Well",
                              positive.ctr, negative.ctr, ...)
        }
        if(!multiplot) {
            if("raw" %in% data.type) {
                if.close.device <- gDevice.new(
                  paste0("boxplot.rawPlate.", object["name"], tag), outpath, ...)
                print(p1)
                gDevice.new(paste0("boxplot.rawWell.", object["name"], tag),
                            outpath, ...)
                print(p2)
            }
            if("norm" %in% data.type) {
                if.close.device <-
                  gDevice.new(paste0("boxplot.normPlate.", object["name"], tag),
                              outpath, ...)
                print(p3)
                gDevice.new(paste0("boxplotnormWell.", object["name"], tag),
                            outpath, ...)
                print(p4)
            }
        } else{
            if(length(data.type) == 2 ) {
                if.close.device <- gDevice.new(
                  paste0("boxplot.raw&norm.", object["name"], tag), outpath, ...)
            } else {
                if.close.device <- gDevice.new(
                  paste0("boxplot.", data.type, ".", object["name"], tag),
                  outpath, ...)
            }
            if(length(data.type) == 2) {
                grid.arrange(p1, p2, p3, p4, nrow = 2 ,ncol = 2)
            }else if (data.type == "raw") {
                grid.arrange(p1, p2, nrow = 2 ,ncol = 1)
            }else if (data.type == "norm") {
                grid.arrange(p3, p4, nrow = 2 ,ncol = 1)
            }
        }
        if(if.close.device) dev.off()
    }
    invisible(list(heatmapRaw = g1, heatmapNorm = g2,
                   boxRawPlate = p1, boxRawWell = p2,
                   boxNormPlate = p3, boxNormWell = p4))
}


## Heatmap of data matrix
## 
## Performs hieratical clustering to data matrix.
## An inner function of \code{\link{cellViz}}
## 
## @param data the data matrix
## @para ... arguments to be passed to \code{cellHeatmap}
## @return Invisibly a list of component by \code{pheatmap}.
## @docType methods
## @import grid
## @import pheatmap
cellHeatmap <- function(data, ...) {
    user.custom <- names(list(...))
    if(!"fontsize" %in% user.custom)
      fontsize <- getOption("opm.heatmap.fontsize")
    if(!"fontsize_row" %in% user.custom)
      fontsize_row <- getOption("opm.heatmap.fontsize_row")
    if(!"fontsize_col" %in% user.custom)
      fontsize_col <- getOption("opm.heatmap.fontsize_col")
    if(!"color" %in% user.custom)
      color <- getOption("opm.heatmap.color")
    if(!"remove.noise" %in% user.custom)
      remove.noise <- getOption("opm.heatmap.rmnoise")

    setHook( "grid.newpage", function()
            pushViewport(viewport(x=1, y=1, width=0.9, height=0.9, name="vp",
                                  just = c("right", "top"))), action="replace" )
    if(remove.noise) {
        medianData <- median(as.vector(data), na.rm = TRUE)
        madData <- mad(as.vector(data), na.rm = TRUE)
        breaks <- seq(medianData - 3 * madData,
                      medianData + 3 * madData, length.out = 51)
    }else{
        breaks <- NA
    }
    g <- pheatmap(data, color = color, fontsize = fontsize,
                  fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                  breaks = breaks, ...)
    setHook("grid.newpage", NULL, "replace")
    grid.text("Well ID", x= -0.02, rot = 90, gp = gpar(fontsize = fontsize))
    grid.text("Plate ID", y= -0.01, gp = gpar(fontsize = fontsize))
    invisible(g)
}

## The boxplot of data of plates and well positions
## An inner function of \code{\link{cellViz}}.
## @param data the data matrix
## @type c("Well", "Plate")
## @positive.ctr matrix
## @negative.ctr matrix
## @... other arguments
## @return Invisibly the ggplot2 function
## @docType methods
## @import ggplot2
cellBoxplot <- function(data, type = c("Well", "Plate"),
                        positive.ctr = NULL, negative.ctr = NULL, ...) {
    if(type == "Plate") {
        data <- suppressMessages(melt(as.data.frame(data)))
    } else if(type == "Well") {
        data <- suppressMessages(melt(as.data.frame(t(data))))
    }
    colnames(data) <- c(type,"Values")
    g <- ggplot(data, aes_string(x = type, y = "Values")) + geom_boxplot()
    if (! is.null(positive.ctr)) {
        if(type == "Plate") {
            data.pc <- suppressMessages(melt(as.data.frame(positive.ctr)))
            colnames(data.pc) <- c(type,"Values")
            g <- g + geom_point(data = data.pc, aes_string(x = type, y = "Values"),
                                size = 1.5, colour = "red")
        } else if(type == "Well") {
            data.pc <- suppressMessages(melt(as.data.frame(t(positive.ctr))))
            g <- g + geom_boxplot(data.pc, aes_string(x = type, y = "Values"),
                                  aes(fill = "red"))
        }
    }
    if (! is.null(negative.ctr)) {
        if(type == "Plate") {
            data.nc <- suppressMessages(melt(as.data.frame(negative.ctr)))
            colnames(data.nc) <- c(type,"Values")
            g <- g + geom_point(data = data.nc, aes_string(x = type, y = "Values"),
                                size = 1.5, colour = "blue")
        } else if(type == "Well") {
            data.nc <- suppressMessages(melt(as.data.frame(t(negative.ctr))))
            g <- g + geom_boxplot(data.nc, aes_string(x = type, y = "Values"),
                                  aes(fill = "blue"))
        }
    }

    invisible(g)
}

#' Hits volcano plot
#' 
#' Visualizes hits by volcano plot.
#'
#' Users can highlight a certain samples during plotting.
#'
#' @param object a cellData object
#' @param outpath diretory of the output figures
#' @param color.highlight a character specifying the color of the hits
#' @param color.background a character specifying the color of the other samples
#' @param highlight.label a vector of characters specifying the names of the
#' samples to be highlighted, with the names are the "barcode:wellID".
#' @param highlight.label.color a character specifying the color of the labels
#' @param ... arguments of the graphic device and ggplot2
#' @docType methods
#' @return  Invisibly an object of \code{ggplot}
#' @examples
#' data(demoCell)
#' op <- options("device")
#' options("device" = "png")
#' labels <- c("Axin1")
#' names(labels) <- c("DSIMGA04:C07")
#' cellSigPlot(oneCell, highlight.label = labels, outpath = tempdir())
#' options(op)
#' @import ggplot2
#' @export
#' 
cellSigPlot <- function(object, outpath = getOption("opm.outpath"),
                        color.highlight = getOption("opm.sig.color.highlight"),
                        color.background = getOption("opm.sig.color.background"),
                        highlight.label = NULL,
                        highlight.label.color = getOption("opm.sig.label.color"),
                        ...) {
    
    if(is.null(object["Sig"]$pvalue)){
        stop("This method requires three or more replicates.")
    }

    objSig <- object["Sig"]
    vc.threshold <- ifelse(apply(objSig$SigMat, 1, any),
                           color.highlight, color.background)
    vc.threshold[is.na(vc.threshold)] <- color.background
    volcanoData <- data.frame(exp = log2(object["qc.data"][, "mean"]),
                              pvalue = -log10(objSig$pvalue),
                              threshold = vc.threshold)

    if.close.device <- gDevice.new(paste(object["name"], "vocalnoPlot",sep="."), outpath, ...)
    exp <- NULL; pvalue <- NULL; label <- NULL
    ## Just for avoiding warnings in Rcheck
    g <- ggplot(data = volcanoData, aes(x = exp, y = pvalue)) +
      geom_point(..., size=2.5, shape=21, colour="black",
                 fill = volcanoData$threshold, alpha = 1) +
                   theme(legend.position = "none") +
                     xlab("log2( intensity )") +
                       ylab("-log10 p-value")
    
    if(!is.null(highlight.label)) {
        TextFrame <- volcanoData[
                          match(names(highlight.label), rownames(volcanoData))
                                 , ]
        if(is.null(TextFrame)) {
            warning("Cannot find the highlight well IDs.
Parameter 'highlight.label' is a vector of labels with the vector names are in
the format of plateID:wellID.")
        } else {
            TextFrame$label <- highlight.label
            TextFrame <- transform(TextFrame,
                                   w = sapply(label, strwidth, "inches"),
                                   h = sapply(label, strheight, "inches"))
            w <- NULL;  h <- NULL      ## Just for avoiding warnings in Rcheck
            g <- g +
              geom_point(..., data = TextFrame, aes(x = exp,y = pvalue),
                         colour = "black", fill = highlight.label.color, 
                         alpha = 1, shape=21, size=3) +
                geom_rect(data = TextFrame,
                           aes(xmin = exp - w/2, xmax = exp + w/2,
                               ymin = pvalue + h/2, ymax = pvalue + 3 * h/2),
                           color = highlight.label.color,
                           fill = highlight.label.color) +
                   geom_text(data = TextFrame,
                             aes(x = exp, y = pvalue, label = label),
                             colour = "black", size = 4, vjust = -1)
        }
    }
    suppressWarnings(print(g))
    if(if.close.device) {
        dev.off()
    }
    invisible(g)
}


#' Hits function analysis
#'
#' Performs function analysis using gProfileR
#'
#' @param object a cellData object
#' @param genemap a data frame, the well-gene specification table
#' @param organism organism name. 
#' @details genemap must include colnames "Barcode","Well","GeneSymbol".
#' organism name can be referred to g:Profiler tool.
#' For example, human: hsapiens, mouse: mmusculus.
#' @param type include both high and low expressed hits or one of them.
#' @param file the filename of the enrichment table (default: disabled)
#' @param ... the arguments of \code{gprofiler}.
#' @return a data frame of the functional report from gProfiler
#' @docType methods
#' @examples
#' data(demoCell)
#' genemap <- read.csv(file.path(system.file("demoData", package = "OperaMate"),
#' "genemap.csv"), stringsAsFactors = FALSE)
#' chart <- cellSigAnalysis(oneCell, genemap, organism = "mmusculus")
#' head(chart)
#' @import gProfileR
#' @export
#'
cellSigAnalysis <- function(object, genemap, organism, type = c("High", "Low"),
                            file = NULL, ...) {

    if (is.null(genemap)) {
        warning("No well-gene specification file provided!")
        return(NULL)
    }
    if(any( !c("Barcode","Well","GeneSymbol") %in% colnames(genemap))){
        warning("Column Barcode, Well and GeneSymbol are not included
                 in the well-gene specification file.")
        return(NULL)
    }
    rownames(genemap) <- paste(genemap$Barcode, genemap$Well, sep = ":")
    objSig <- object["Sig"]$SigMat
    sig <- objSig[match(rownames(genemap), rownames(objSig)), ]
    ID <- NULL
    if ("High" %in% type) 
      ID <- c(ID, rownames(objSig)[objSig[, "High"]])
    if ("Low" %in% type) 
      ID <- c(ID, rownames(objSig)[objSig[, "Low"]])
    query <- genemap[match(na.omit(ID), rownames(genemap)), "GeneSymbol"]
    pars <- list(...)
    gprofile.pars <- formals(gprofiler)
    pars <- pars[names(pars) %in% names(gprofile.pars)]
    gprofile.chart <- try(gprofiler(query, organism = organism, pars))
    if (class(gprofile.chart) == "try-error") {
        warning("Function gProfile failed.")
        return(NULL)
    }

    if (!is.null(file)) {
        write.table(gprofile.chart, file = file,
                    row.names = FALSE, col.names = TRUE,
                    quote = FALSE, sep = "\t")
    }
    return(gprofile.chart)
}

