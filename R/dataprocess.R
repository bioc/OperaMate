#' Data process and analysis pipeline
#' 
#' A systematical pipeline for opera data importing, normalization, quality
#' control, signature gene detection and visualization.
#' 
#' @param cellformat character specifying the format of the reports
#' @param datapath character specifying the location of the reports
#' @param platemap data frame, required when cellformat = "Others".
#' See an example as \code{\link{platemap}}.
#' @param Prefix.Barcode character, the Barcode prefix
#' @param genemap data frame, the mapping file of all samples.
#' See an example as \code{\link{genemap}}.
#' @param cellnames vector of character specifying the data types.
#' The types are the parameters in the Columbus system reports.
#' @param ctrwell a vector of characters specifying the control well IDs
#' @param expwell a vector of characters specifying the control sample IDs
#' @param norm.method norm.method = c("Both","Plate","Z","Ctr","None")
#' @param gDevice the graphic device (default: png)
#' @param outpath a character string naming the location the figures to
#' generate to
#' @param BatcheffectSample c("all", "exp", "ctr"), "exp" and "ctr" consider
#' only the sample and control wells respectively.
#' @param BatcheffectType  c("heatmap","median","PCA-plate","PCA-well")
#' @param BatcheffectifCompare logic, comparison before and after normalization
#' will be shown if the data normalization has been processed.
#' @param qth numeric, the threshold in quality control
#' @param Sig.method c("ktsd","ksd","kmsd")
#' @param Sig.threshold numeric, the threshold for hit identification
#' @param anno.file the path to generate the report
#' @param viz.par list, arguments to be passed to \code{\link{cellViz}}
#' @param email the email for the DAVID webservice registration
#' @param david.terms a list of DAVID annotation categories
#' @param david.count numeric, number of items
#' @param david.threshold numeric, threshold in the DAVID analysis
#' @param david.termN numeric, number of terms to show in the figure
#' @param david.type c("all","high","low")
#' @param threshold.pvalue numeric, threshold of pvalues in the t-test between
#' the sample and control replicates
#' @param ... arguments of the graphic device
#' @return a list of two components:  a list of cellData objects and  the
#' annotated table of each well
#' @export
#' @examples
#' data(platemap); data(genemap)
#' platemap$path <- file.path(system.file("Test",package = "OperaMate"),
#'                            platemap$path)
#' (outpath <- tempdir())
#' operaReport <- operaMate(cellformat = "Others",
#' platemap = platemap,genemap = genemap,
#' Prefix.Barcode = "DSIMGA",outpath=outpath,
#' expwell = paste0(rep(LETTERS[2:15],each=20),
#' rep(formatC(3:22,width=2,flag=0),times=14)),
#' ctrwell = c("C23","E23","G23"),
#' email = "cliu@@sjtu.edu.cn",
#' norm.method="Both",BatcheffectSample="exp",
#' Sig.method="ktsd",david.type="all",
#' width=960,height=960)
#' head(operaReport$data.anno)
#' 
operaMate <- function(cellformat=c("Tab","Matrix","Others"),datapath="./",
                      platemap=NULL, genemap=NULL,
                      cellnames=c("Average.Intensity.of.Nuclei",
                        "Average.Intensity.of.Cytoplasm",
                        "Average.Total.Intensity","Average.Intensity.Ratio"),
                      ctrwell = NULL,
                      expwell = NULL,
                      Prefix.Barcode = "DSIMGA",
                      norm.method=c("Both","Plate","Z","Ctr","None"),
                      gDevice="png",
                      outpath="./",
                      BatcheffectSample = c("all", "exp", "ctr"),
                      BatcheffectType=c("heatmap","median","PCA"),
                      BatcheffectifCompare=TRUE,
                      viz.par = list(),
                      qth = 0.05,
                      Sig.method=c("ktsd","ksd","kmsd"),
                      Sig.threshold = NULL,
                      threshold.pvalue=0.05,
                      anno.file = file.path(outpath,"OperaMateReport.txt"),
                      email = NA,
                      david.type = c("all","high","low"),
                      david.terms = NULL,
                      david.count=2L,
                      david.threshold=0.1,
                      david.termN=NA,
                      ...
                      ){
    if(!file.exists(outpath)){
        dir.create(outpath)
    }
    message("[",format(Sys.time(), "%m-%d-%Y %T"),"]")
    message(" OperaMate Data Processing & Analysis")
    message("********************************************************")
    message("Loading data ...")
    lstPlates <- loadAll(cellformat=cellformat,datapath=datapath,
                         platemap=platemap,prefix=Prefix.Barcode)
    cellnames.1 <- cellnames[cellnames %in% names(lstPlates[[1]]["data"])]
    if("Average.Total.Intensity" %in% cellnames){
        cellnames.1 <- c("Average.Intensity.of.Nuclei",
                         "Average.Intensity.of.Cytoplasm", cellnames.1)
    }
    lstCells <- list()
    for(cellname in cellnames.1){
        oneCell <- cellData(cellname,ctrwell=ctrwell,expwell=expwell)
        lstCells[[cellname]] <- cellLoad(oneCell,lstPlates)
    }

    message("Data normalization ...")
    lstCells <- lapply(lstCells,function(cell){
        cellNorm(cell,norm.method=norm.method)
    })
    message("Data biases visualization ...")
    for(cell in lstCells){
        for(type in BatcheffectType){
            cellViz(cell,exps=BatcheffectSample,type=type,gDevice=gDevice,
                    outpath=outpath,ifCompare=BatcheffectifCompare,
                    control.par=viz.par,...)
        }
    }
    message("Quality control ...")
    lstCells <- lapply(lstCells,function(cell){
        cellQC(cell, qth=qth,
               outpath=outpath,
               gDevice=gDevice,...)
    })
    if("Average.Total.Intensity" %in% cellnames){
        lstCells[["Average.Total.Intensity"]] <-
          cellMean(lstCells[["Average.Intensity.of.Nuclei"]],
                   lstCells[[ "Average.Intensity.of.Cytoplasm"]],
                   "Average.Total.Intensity")
    }
    lstCells <- lstCells[cellnames]
    message("Hit detection ...")
    lstCells <- lapply(lstCells,function(cell){
        cell <- cellSig(cell, method=Sig.method,
                        threshold.method = Sig.threshold,
                        threshold.pvalue = threshold.pvalue,
                        adjust.method = "fdr",
                        outpath=outpath,
                        gDevice=gDevice,...)
        cellSigplot(cell,
                    outpath=outpath,
                    gDevice=gDevice,...)
        return(cell)
    })
    message("Annotation ... ")
    data.anno <- GenerateReport(lstCells,genemap,file=anno.file,
                                outpath=outpath,
                                gDevice=gDevice,...)
    message("Hit analysis ...")
    david.anno <- lapply(lstCells,function(cell){
        lapply(david.type,function(type){
            cellSigAnalysis(cell,genemap,email,
                            type=type,
                            david.terms = david.terms,
                            count=david.count,
                            threshold=david.threshold,
                            showTermNum=david.termN,
                            file=file.path(outpath,
                              paste0("DAVID-",cell["name"],"-",type,".txt")),
                            outpath=outpath,gDevice=gDevice,...)
        })
    })

    return(list(lstCells = lstCells,data.anno=data.anno))
}
