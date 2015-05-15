## cellData: data of each type.
##
#' The cellData class
#'
#' The main class used in \code{OperaMate} to hold all levels of experiment
#' data of a specific type.
#' 
#' @slot name character, one parameter in the Columbus system report.
#' @slot ctrwell a character vector, the control well IDs, e.g. B05.
#' @slot expwell a character vector, the sample well IDs, e.g. C15.
#' @slot origin.data a numeric matrix,
#' the raw data matrix with rows the well IDs and columns the plate IDs.
#' @slot norm.data a numeric matrix, the normalized data.
#' @slot qc.data a numeric matrix, the data after quality
#' control, with the rows are "barcode:wellID" and columns are the data of
#' all replicated samples and their means, the p-values in quality control
#' and if they have passed the quality control.
#' @slot norm.method character the normalization method.
#' @slot qc.th numeric, the threshold in the quality control.
#' @slot plate.quality a numeric matrix, the quality data with
#' the rows are the barcode and columns are the replicateIDs.
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
#' \code{cellData(name, ctrwell, expwell = NULL)}. \code{name}: character,
#' one parameters in the Columbus system reports.
#' }
#' \item{Show}{signature(object="cellLoad"). Displays object content as text.}
#' \item{Accessor}{x[i]. \code{x}: a cellData object;
#' \code{i}: character, a cellData slot name.}
#' }
#' @examples
#' oneCell <- cellData(name="Average Intensity of Nuclei",
#' ctrwell = c("C23","E23","G23"))
#' oneCell
#' oneCell["name"]
#' @exportClass cellData
#' @docType class
#' @name cellData-class
#' @rdname cellData-class
#' @aliases cellData,cellData-class
setClass("cellData", slots = c(
                       name    = "character",
                       ctrwell = "character",
                       expwell = "character",
                       origin.data = "matrix",
                       norm.data = "matrix",
                       qc.data = "data.frame",
                       norm.method = "character",
                       qc.th = "numeric",
                       plate.quality = "matrix",
                       Sig = "list"
                       ),
         prototype=list(
           qc.th = 0.05,
           norm.method = c("Both","Plate","Well","Ctr","Z","None")
           )
         )
#' @name cellData
#' @rdname cellData-class
#' @aliases cellData,cellData-method
#' @export
cellData <- function(name,ctrwell,expwell=NULL,...){
    if(length(intersect(ctrwell,expwell))){
        stop("Intersection of control and sample well IDs exist")
    }
    name <- gsub(" ",".",name)
    args <- list(...)
    if(is.null(args[["norm.method"]])){
        norm.method <- c("Both","Plate","Well","Ctr","Z","None")
    }else if(!args[["norm.method"]]%in%
             c("Both","Plate","Well","Ctr","Z","None")){
        stop("norm.method must be one of (Both,Plate,Well,Ctr,Z,None)")
    }else{
        norm.method <- args[["norm.method"]]
    }
    if(is.null(args[["qc.th"]])){
        qc.th <- 0.05
    }else{
        qc.th <- args[["qc.th"]]
    }
    if(is.null(expwell)){
        expwell <- character(0)
    }
    new("cellData",name=name,ctrwell=ctrwell,expwell=expwell,
        norm.method=norm.method,qc.th=qc.th)
}

## Getter
#' @rdname cellData-class
setMethod("[", signature=c(x="cellData",i="character"),function(x,i){
    if(i=="name"){return(x@name)}else{}
    if(i=="origin.data"){return(x@origin.data)}else{}
    if(i=="norm.data"){return(x@norm.data)}else{}
    if(i=="qc.data"){return(x@qc.data)}else{}
    if(i=="plate.quality"){return(x@plate.quality)}else{}
    if(i=="Sig"){return(x@Sig)}else{}
    if(i=="norm.method"){return(x@norm.method)}else{}
    if(i=="ctrwell"){return(x@ctrwell)}else{}
    if(i=="expwell"){return(x@expwell)}else{}
})
## Show method
#' @rdname cellData-class
#' @aliases show,cellData-method
setMethod("show",signature=(object="cellData"),function(object){
    cat("An object of cellData class.\n")
    cat("Parameter: ", object@name,".\n")
    cat("Raw Data: ")
    if(length(object@origin.data)){
        str(object@origin.data)
        cat("\n")
    }else{
        cat("NULL.\n")
    }
    cat("Data Normalization: ",
        ifelse(length(object@norm.data),"Done","ToDo..."),";\n",sep="")
    cat("     method: ",
        paste(object@norm.method, collapse=" "),
        ".\n",sep="")
    cat("Quality Control: ",
        ifelse(length(object@qc.data),"Done","ToDo..."),".\n",sep="")
    cat("Significant hits: \n")
    if(length(object@Sig)){
        show(object@Sig$stats)
    }else{
        cat("      ToDo...\n")
    }
})

#' Data importing
#' 
#' Extracts data of a specific type in a list of \code{expData} objects to
#' initialize a \code{cellData} object.
#' @param object a cellData object
#' @param lstPlates a list of expData objects
#' @return a \code{cellData} object, with initialized slot \code{origin.data}
#' @docType methods
#' @examples
#' data(platemap)
#' platemap$path <- file.path(system.file("Test",package = "OperaMate"),
#' platemap$path)
#' lstPlates <- loadAll(cellformat="Tab",datapath=dirname(platemap$path[1]))
#' data(demoCell)
#' oneCell <- cellLoad(oneCell, lstPlates)
#' str(oneCell["origin.data"])
#' @export
#' @rdname cellLoad
setGeneric("cellLoad",function(object,lstPlates){
    standardGeneric("cellLoad")
})
#' @rdname cellLoad
#' @aliases cellLoad,cellData-method
#' 
setMethod("cellLoad", signature = "cellData", function(object,lstPlates){
    allWell <- c( object@expwell, object@ctrwell )
    if( is.null(object@expwell) ){
        allWell <- unique(unlist(lapply(lstPlates,
                                        function(onePlate) onePlate["wellID"])))
        object@expwell <- allWell[!allWell %in% object@ctrwell]
    }
    onecelldata <- matrix(NA,length(allWell),length(lstPlates))
    i <- 0
    for(onePlate in lstPlates){
        i <- i + 1
        tmpdata <- (onePlate["data"])[[object@name]]
        onecelldata[,i] <- tmpdata[match(allWell,
                                         onePlate["wellID"])]
    }
    colnames(onecelldata) <- sapply(lstPlates,function(x) x["name"])
    rownames(onecelldata) <- allWell
    object@origin.data <- onecelldata
    return(object)
})

#' Data normalization
#' 
#' Normalizes raw data based on different normalization methods.
#' 
#' Method description: "Both" employes the median polish algorithm which
#' divides data by the median of their plates and wells recursively, while
#' "Plate" only divides data by the median of their plates;"Z" substracts
#' data by their plate medians, and then divides by the median absolute
#' deviations; "Ctr" divides data by the mean of their plate controls;
#' "None" avoids the data normalization in this step. The first three
#' methods are based on the assumption that most samples display no
#' biological effects in the assay be analyzed. They are often more
#' effective than "Ctr" method as to the high throughput screening.
#'
#' @param object a cellData object
#' @param norm.method norm.method=c("Both","Plate","None","Z","Ctr")
#' 
#' @return a \code{celldata} object with initialized slot \code{norm.data}
#' @examples
#' data(demoCell)
#' oneCell <- cellNorm(oneCell, norm.method="Plate")
#' str(oneCell["norm.data"])
#' @export
#' @rdname cellNorm

setGeneric("cellNorm",function(object,norm.method){
    standardGeneric("cellNorm")
})
#' @rdname cellNorm
#' @aliases cellNorm,cellData-method
setMethod("cellNorm", signature = "cellData",
          function(object,norm.method=c("Both","Plate","None","Z","Ctr")){
    if(length(object@norm.method)>1){
        object@norm.method <- norm.method[1]
    }else if(length(norm.method)==1&&object@norm.method!=norm.method){
        object@norm.method <- norm.method
    }else{}
    norm.method <- object@norm.method
    if( norm.method == "None" ){
        object@norm.data <- object@origin.data
    }else if( norm.method == "Plate" ){
        object@norm.data <- t( t(object@origin.data)
                              /apply(object@origin.data,2,median,na.rm=TRUE) )
    }else if( norm.method == "Both" ){
        log <- capture.output(
          {tmp <- medpolish(log2(object@origin.data),na.rm=TRUE)})
        object@norm.data <- 2^(tmp$residuals)
    }else if( norm.method == "Ctr" ){
        if(is.vector(object@ctrwell)){
            lstctr <- apply( object@origin.data, 2, function(x)
                            mean(x[object@ctrwell],na.rm=TRUE) )
        }else if (is.list(object@ctrwell)){
            lstctr <- apply(names(object@ctrwell),function(ctr){
                mean(object@origin.data[object@ctrwell[[ctr]],ctr],na.rm=TRUE)})
        }
        object@norm.data <- t( t(object@origin.data)/lstctr )
    }else if( norm.method == "Z"){
        object@norm.data <- t( (t(object@origin.data)-
                                apply(object@origin.data,2,median,na.rm=TRUE))
                              /apply(object@origin.data,2,mad,na.rm=TRUE) )
    }else{
        stop("Set your own normalization rule to object@norm.method
             and the corresponding normalized data to object@norm.data.")
    }
    return( object )
})


#' Quality control
#' 
#' Checks quality of all plates and then wells.
#'
#' Requires three or more replicated samples.
#' 
#' @param object a cellData object
#' @param qth a numeric threshold of pvalues.
#' @param outpath a character string naming the location the figure to
#' generate to.
#' @param gDevice the graphic device (default: dev.new).
#' gDevice=NA generates no figure.
#' @param ... arguments of the graphic device
#' 
#' @return a \code{cellData} object with intialized slot \code{qc.data}
#' and \code{plate.quality}
#' @rdname cellQC
#' @examples
#' data(demoCell)
#' oneCell <- cellQC( oneCell, qth = .05, gDevice = NA)
#' str(oneCell["qc.data"])
#' str(oneCell["plate.quality"])
#' @export
#' 
setGeneric("cellQC",function(object,
                             qth=0.05, outpath="./", gDevice="dev.new",...){
    standardGeneric("cellQC")
})
#' @rdname cellQC
#' @aliases cellQC,cellData-method
setMethod("cellQC", signature = "cellData",
          function(object, qth=0.05, outpath="./", gDevice="dev.new",...){
    if(object@qc.th!=qth){
        object@qc.th <- qth
    }
    out <- plateQC(object@norm.data,qth=qth)
    plate.quality <- out$plate.quality
    Qplate.intensity <- out$Qplate.intensity
    object@plate.quality <- plate.quality
    if(ncol(plate.quality)==1){
        colnames(Qplate.intensity) <- sapply(
          strsplit(colnames(Qplate.intensity),"-"),function(x) x[[1]])
        gene.quality <- as.vector(Qplate.intensity)
        names(gene.quality) <- paste(rep(colnames(Qplate.intensity),
                                         each=nrow(Qplate.intensity)),
                                     rownames(Qplate.intensity),sep="-")
        object@qc.data <- as.data.frame(gene.quality)
    }else{
        gene.quality <- wellQC(Qplate.intensity,qth=qth,
                               prefix = object@name,outpath=outpath,
                               gDevice=gDevice,...)
        object@qc.data <- gene.quality
    }
    return(object)
})

#' Hit identification
#' 
#' Detects samples those are most different from the negative controls.
#' 
#' @param object a cellData object
#' @param method method=c("ktsd","ksd","kmsd").
#' Details are referred in the vignette.
#' @param threshold.method numeric
#' @param digits integer, the number of digits used to show the results
#' @param outpath a character string naming the location the figure to
#' generate to. Only for method = "ktsd".
#' @param gDevice the graphic device (default: dev.new),
#' gDevice=NA generates no figure.  Only for method = "ktsd".
#' @param threshold.pvalue numeric, threshold of pvalues in the t-test
#' between the sample and control replicates
#' @param adjust.method pvalue correction method
#' @param ... arguments of the graphic device
#' @return a \code{cellData} object with initialized slot \code{Sig}.
#' @rdname cellSig
#' @examples
#' data(demoCell)
#' oneCell <- cellSig( oneCell, method = "ktsd",gDevice=NA )
#' names(oneCell["Sig"])
#' 
#' @export
#' 
setGeneric("cellSig",function(object,method=c("ktsd","ksd","kmsd"),
                              threshold.method=NULL,digits=3,outpath="./",
                              gDevice="dev.new",threshold.pvalue=0.05,
                              adjust.method=p.adjust.methods,... ){
    standardGeneric("cellSig")
})
#' @rdname cellSig
#' @aliases cellSig,cellData-method
setMethod("cellSig", signature = "cellData",
          function(object, method=c("ktsd","ksd","kmsd"),
                    threshold.method=NULL,digits=3,outpath="./",
                    gDevice="dev.new",threshold.pvalue=0.05,
                    adjust.method=p.adjust.methods,... ){
    flag <- 1
    if(is.null(threshold.method)){
        if(method == "ktsd"){
            threshold.method <- 0.05
        }else if(method =="ksd"||method == "kmsd"){
            threshold.method <- 3
        }
    }
    goodWellID <- which(object@qc.data[,"pval.pass"])
    repdata <- object@qc.data[ goodWellID,
                              1:(ncol(object@qc.data)-3) ]
    if(is.null(ncol(repdata)) || ncol(repdata)<3){
        message("T-test method will not be performed ")
        message("since it requires three or more replicate data.")
        flag <- 0
        t.pvalue <- NULL
    }else{
        sig0 <- sig.t(repdata,object@ctrwell,
                      pth=threshold.pvalue,adjust.method=adjust.method)
        t.pvalue <- sig0$pvalue[match(rownames(object@qc.data),
                                      rownames(sig0$sig.matrix))]
        names(t.pvalue) <- rownames(object@qc.data)
    }
    gene.exp <- object@qc.data[ goodWellID
                               ,"mean" ]
    names(gene.exp) <- rownames(object@qc.data[goodWellID,])
    if(method == "ktsd"){
        sig <- sig.ktsd(gene.exp, 
                        pth = threshold.method,
                        prefix = object@name, gDevice=gDevice,
                        outpath=outpath,...)
    }else if (method == "ksd"){
        sig <- sig.ksd(gene.exp,k=threshold.method, digits=digits)
    }else if (method == "kmsd"){
        sig <- sig.kmsd(gene.exp,k = threshold.method, digits=digits)
    }
    threshold <- c(threshold=threshold.method,sig$threshold)
    if(flag){
        SigMat <- sig$sig.matrix & sig0$sig.matrix
    }else{
        SigMat <- sig$sig.matrix
    }
    stats <- apply(SigMat,2,sum,na.rm=TRUE)
    SigMat <- SigMat[match(rownames(object@qc.data),
                           rownames(SigMat)),]
    rownames(SigMat) <- rownames(object@qc.data)
    object@Sig <- list(SigMat = SigMat,
                       threshold = threshold,
                       stats = stats,
                       pvalue = t.pvalue
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
#' meanCell <- cellMean(oneCell,oneCell,"meanCell")
#' meanCell
setGeneric("cellMean", function(cell1,cell2,name){
    standardGeneric("cellMean")
})
#' @rdname cellMean
#' @aliases cellMean,cellData-method
setMethod("cellMean", signature = c(cell1="cellData", cell2="cellData",name="character"),
          function(cell1,cell2,name){
              oneCell <- cellData(name = name,
                             expwell = cell1@expwell,
                             ctrwell = cell1@ctrwell,
                             qc.th = unique(c(cell1@qc.th,cell2@qc.th)),
                             norm.method =
                             unique(c(cell1@norm.method,cell2@norm.method)))
              oneCell@origin.data <- 
                (cell1@origin.data + cell2@origin.data)/2
              oneCell@norm.data <-
                (cell1@norm.data + cell2@norm.data )/2
              oneCell@plate.quality <- 
                (cell1@plate.quality & cell2@plate.quality)
              oneCell@qc.data <- data.frame(
                (cell1@qc.data[,1:(ncol(cell1@qc.data)-2)] 
                 + cell2@qc.data[,1:(ncol(cell2@qc.data)-2)])/2,
                pval = ifelse(
                  cell1@qc.data[,"pval"]>cell2@qc.data[,"pval"],
                  cell1@qc.data[,"pval"],
                  cell2@qc.data[,"pval"]),
                pval.pass = cell1@qc.data[,"pval.pass"]&
                cell2@qc.data[,"pval.pass"])
              return(oneCell)
          }
          )
