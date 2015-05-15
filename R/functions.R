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
#' \item{\code{file.name}:}{character, the name of the report.}
#' \item{\code{format.type}:}{character, only ''Tab'' and ''Matrix'' are supported in
#' the current version.}
#' \item{\code{Barcode}:}{character, the barcode of the plates.}
#' \item{\code{rep.id}:}{character, the ID to distinguish the replicated plates.}
#' \item{\code{path}:}{character, the full path of the report.}
#' }
#' }
#' }
#' @examples
#' data(platemap)
#' str(platemap)
#' @return platemap: a data frame
NULL
#' @name genemap
#' @aliases genemap
#' @rdname data
#' @section genemap:
#' \describe{
#' \item{Description}{The gene names and informations of each well.
#' The correponding file is generated during the design of assays.
#' Required for function \code{\link{GenerateReport}}
#' and \code{\link{cellSigAnalysis}}.}
#' \item{Format}{data.frame with the following required column names:
#' \describe{
#' \item{\code{Barcode}:}{character, the barcode of the plates.}
#' \item{\code{Well}:}{character, the well ID.}
#' \item{\code{GeneSymbol}:}{character, the annotated gene names of the well.}
#' \item{\code{Gene.ID}:}{character, Entrez gene ID.}
#' }}
#' }
#' @examples
#' data(genemap)
#' str(genemap)
#' @return genemap:a data frame 
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

## Quality control of each sample
##
## Performs quality control to each well based on its duplicated samples.
## An internal function of \code{\link{cellQC}}.
## @param intensity matrix, with rows the wells and columns the plates
## @param qth numeric, the theshold of the pvalues
## @param prefix character, the name of generated figure
## @param outpath a character string naming the location the figures to
## generate to
## @param gDevice the graphic device (default: dev.new)
## @param ... arguments of the graphic device
## @return a matrix containing the data after quality control, with the rows are
## "barcode:wellID" and columns are the data of each replicated samples and
## their means, the p-values in quality control and if they have passed the
## quality control.
## 
wellQC <- function(intensity, qth=0.05,
                   prefix="", outpath="./", 
                   gDevice="dev.new",...){
    splitcol <- lapply(split(colnames(intensity),
                             unlist(lapply(colnames(intensity),function(x)  {
                                 unlist(strsplit(x,"-"))[1]}))),sort)
    gene.intensity <- do.call(rbind,lapply(names(splitcol),function(x){
        sub.data <- intensity[,splitcol[[x]]]
        rownames(sub.data) <- paste(x,rownames(sub.data),sep=":")
        return(sub.data)
    }))
    colnames(gene.intensity) <- paste("s",c(1:ncol(gene.intensity)),sep="")
    mean.intensity <- apply(gene.intensity,1,mean,na.rm=TRUE)
    sd.intensity <- apply(gene.intensity,1,sd,na.rm=TRUE)
    mean.transform <- ifelse(mean.intensity>1,mean.intensity,1/mean.intensity)
    mean.th <- quantile(mean.transform,0.99,na.rm=TRUE)
    useID <- (!is.na(sd.intensity) & mean.transform<mean.th)
    fit.mean <- mean.transform[useID]
    fit.sd <- sd.intensity[useID]
    delID2 <- is.na(mean.transform)
    fit.intensity <- smooth.spline(fit.mean,fit.sd)
    fit.dist <- (sd.intensity[!delID2]-
                 predict(fit.intensity,mean.transform[!delID2])$y)/
                   mean.transform[!delID2]
    fit.dist.1 <- c(fit.dist[fit.dist<0],(-1)*fit.dist[fit.dist<0])
    pval <- 1-pnorm(fit.dist,0,sd(fit.dist.1,na.rm=TRUE))
    pval.f <- rep(NA,length(mean.intensity))
    pval.f[!delID2] <- pval; pval <- pval.f
    pval.pass <- (pval > qth)
    gene.quality <- data.frame(gene.intensity, mean=mean.intensity,
                               pval, pval.pass, stringsAsFactors=FALSE)
    if(!is.na(gDevice)){
        if(gDevice=="dev.new"){
            dev.new();
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice)
                 )(
                   file.path(outpath,
                             paste0(prefix, ".replicates.sd2mean.",gDevice)),...
                   )
        }
        par(mfrow=c(2,1))
        plot(smooth.spline(fit.sd,fit.mean),xlab="SD",ylab="Mean",
             main="Curve of replicates SD to mean",type="b",lwd=2)
        plot(density(fit.dist[!is.na(fit.dist)]),
             main="Replicates SD difference between ideal and experiment",
             xlab=paste("Difference, sd = ",
               formatC(sd(fit.dist,na.rm=TRUE),width=2)))
        abline(v=qnorm(1-qth,0,sd(fit.dist.1,na.rm=TRUE)),col="red")
        if(gDevice!="dev.new"){
            graphics.off();
        }
    }
    return(gene.quality)
}

## Quality control of duplicated plates
##
## Performs quality control based on duplicated plates.
## An internal function of \code{\link{cellQC}}.
## @param intensity matrix, with rows the intensity of the wells and columns
## the intensity of the plates
## @param qth numeric, the theshold of pvalue.
## @return a list with two components: an intensity matrix after quality control
## and the quality logical matrix.
plateQC <- function(intensity, qth=0.05){
    splitcol <- split(colnames(intensity),
                      unlist(lapply(colnames(intensity),
                                    function(x){ unlist(strsplit(x,"-"))[1]})))
    plate.test <- do.call(rbind,lapply(splitcol,function(coln){
        if(length(coln)==1){
            return(1)
        }else{
            out <- intensity[,coln]
            comb.col <- combn(1:length(coln),2)
            pval.vec <- apply(comb.col,2,function(colvec){
                return( t.test(out[,colvec[1]],out[,colvec[2]])$p.value)
            })
            return(pval.vec)
        }}))
    if(ncol(plate.test)==1){
        plate.quality <- plate.test
        mode(plate.quality) <- "logical"
        colnames(plate.quality) <- "s1"
        Qplate.intensity <- intensity
    }else{
        Qpass <- apply(plate.test>qth,1,function(
          qwell){
            tmp <- matrix(0,length(qwell),length(qwell))
            tmp[upper.tri(tmp)] <- qwell
            tmp[lower.tri(tmp)] <- qwell
            apply(tmp,1,function(x) any(x==1,na.rm=TRUE))})
        plate.quality <- t(Qpass)
        colnames(plate.quality) <- paste("s",1:ncol(plate.quality),sep="")
        Qplate.intensity <- do.call(cbind,lapply(colnames(Qpass),function(x){
            qID <- Qpass[,x]
            sub.data <- intensity[,splitcol[[x]]]
            if(sum(qID)==length(qID)){
            }else if(sum(qID)>1){
                mean.val <- apply(sub.data[,qID],1,mean,na.rm=TRUE)
                sub.data[,!qID] <- mean.val
            }else if (sum(qID)==1){
                mean.val <- sub.data[,qID]
                sub.data[,!qID] <- mean.val
            }else if (sum(qID)==0){
                sub.data <- matrix(NA,nrow(sub.data),ncol(sub.data))
            }
            colnames(sub.data) <- splitcol[[x]]
            return(sub.data)
        }))
    }
    return(list(Qplate.intensity=Qplate.intensity,
                plate.quality=plate.quality))
}

## Significant samples detection by ksd method
##
## An internal function of \code{\link{cellSig}}.
## @param gene.exp vector, the data of each sample
## @param k numeric, hits are the samples beyond the range
## [mean-k*sd, mean+k*sd].
## @param digits integer, the number of digits used to show the results
## @keywords internal
## @return a list including the hit matrix and the thresholds 
sig.ksd <- function(gene.exp,k=3,digits=3){
    m.exp   <- mean(gene.exp,na.rm=TRUE)
    sd.exp  <- sd(gene.exp,na.rm=TRUE)
    ths  <- c( Low = m.exp - k * sd.exp,
              High = m.exp + k * sd.exp )
    sig.matrix <- cbind( gene.exp < ths["th.low"],
                        gene.exp > ths["th.high"])
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- names(gene.exp)
    colnames(sig.matrix) <- c("Low","High")
    ths <- round(ths,digits=digits)
    return(list(sig.matrix=sig.matrix,threshold=ths))
}
## Significant samples detection by kmsd method
##
## An internal function of \code{\link{cellSig}}.
## @param gene.exp vector, the data of each sample
## @param k numeric, hits are the samples beyond the range
## [median-k*mad, median+k*mad].
## @param digits integer, the number of digits used to show the results
## @return a list including the hit matrix and the thresholds
## 
sig.kmsd <- function(gene.exp,k=3,digits=3){
    m.exp   <- median(gene.exp,na.rm=TRUE)
    sd.exp  <- mad(gene.exp,na.rm=TRUE)
    ths  <- c ( Low = m.exp - k * sd.exp,
               High = m.exp + k * sd.exp )
    sig.matrix <- cbind( gene.exp < ths["th.low"],
                        gene.exp > ths["th.high"])
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- names(gene.exp)
    colnames(sig.matrix) <- c("Low","High")
    ths <- round(ths,digits=digits)
    return(list(sig.matrix=sig.matrix,threshold=ths))
}
## Significant samples detection by ktsd method
##
## An internal function of \code{\link{cellSig}}.
## @param gene.exp vector, the data of each sample
## @param pth numeric, the threshold of t-scores
## @param digits integer, the number of digits used to show the results
## @param outpath a character string naming the location the figure to
## generate to
## @param gDevice the graphic device (default: dev.new)
## @param prefix the title of the figure
## @param ... arguments of the graphic device
## @return a list including the hit matrix and the thresholds
## 
sig.ktsd <- function(gene.exp,
                     pth = 0.05, digits=3, 
                     outpath="./",prefix="",gDevice="dev.new",...){
    fit.exp <- fitT(gene.exp,
                    prefix=prefix,
                    outpath=outpath,
                    gDevice=gDevice,...)
    th <- qt(1-pth/2,df=fit.exp$estimate["df"])
    ths <- c( Low = (-1)*th*fit.exp$estimate["s"] + fit.exp$estimate["m"],
             High = th*fit.exp$estimate["s"] + fit.exp$estimate["m"] )
    sig.matrix <- cbind( gene.exp < ths[1],
                        gene.exp > ths[2])
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- names(gene.exp)
    colnames(sig.matrix) <- c("Low","High")
    ths <- round(ths,digits=digits)
    names(ths) <- c("Low","High")
    return(list(sig.matrix=sig.matrix,threshold=ths))
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
sig.t <- function(repdata, ctrwell,
                  pth = 0.05, adjust.method=p.adjust.methods){
    lst.plate.info <-
      split(repdata, sapply(strsplit(rownames(repdata),":"),function(x) x[1]))
    pval.well <- unlist(lapply(lst.plate.info,function(pmat){
        ctrexp <-
          as.vector( pmat[grep(paste(ctrwell,collapse="|"), rownames(pmat)), ] )
        apply(pmat,1,function(exp) {
            ifelse(any(is.na(exp)),NA,
                   t.test(exp, ctrexp,alternative="less")$p.value) })
    }))
    names(pval.well) <- sapply(strsplit(names(pval.well),"[.]"),function(x)x[2])
    pval.well.twoside <- ifelse(pval.well<0.5, pval.well,
                                1-pval.well)*2
    padjust.well <- p.adjust(pval.well.twoside,adjust.method)
    sig.ID <- (padjust.well < pth)
    sig.matrix <- cbind( pval.well<0.5 & sig.ID,
                        pval.well>0.5 & sig.ID)
    mode(sig.matrix) <- "numeric"
    rownames(sig.matrix) <- rownames(repdata)
    colnames(sig.matrix) <- c("Low","High")
    return(list(sig.matrix=sig.matrix,
                pvalue = padjust.well))
}
## Fit data to student's t distribution
## 
## Fits data to a student's t distribution.
## @param gene.exp vector
## @param plot logic, if to plot the figure
## @param prefix character, the name of the figure
## @param outpath a character string naming the location the figure to
## generate to. 
## @param gDevice the graphic device (default: dev.new),
## gDevice=NA generates no figure.  
## @param ... arguments of the graphic device
## @param xlab character, X axis label
## @param ylab character, Y axis label
## @return  An object of class "fitdistr" referred in MASS package

fitT <- function( gene.exp, plot=TRUE,
                 prefix="", outpath="./",
                 gDevice="dev.new",xlab="intensity",ylab="t distribution",...){
    naID <- is.na(gene.exp)
    gene.exp[naID] <- sample(gene.exp[!naID],sum(naID),replace=TRUE)
    suppressWarnings( fit <- fitdistr( gene.exp, "t" ) )
    if(!is.na(gDevice)){
        if(gDevice=="dev.new"){
            dev.new();
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice)
                 )(file.path(outpath,
                             paste0(prefix, "IntensityDistribution.",gDevice)),
                   ...)
        }
        fit.gene.exp <- (gene.exp-fit$estimate["m"]) / fit$estimate["s"]
        qqplot( fit.gene.exp,
               rt(1000, df = fit$estimate["df"]),
               xlab = xlab,
               ylab = ylab,
               main = "QQ plot with t distribution")
        abline(c(0,1))
        if(gDevice!="dev.new"){
            graphics.off();
        }
    }
    return(fit)
}

## Visualize data of one plate
##
## An internal function of \code{\link{plateViz}}.
## 
## @param onePlateData a vector containing the data of each well and
## the names are the well IDs
## @param main character, the title of the heatmap
## @param ... auguments to be passed to \code{pheatmap}
## @return Invisibly a list of component by \code{pheatmap}.
## 
## @import pheatmap
## 
showPlate <- function(onePlateData,main=NULL, ...){
    tmp <-  split(onePlateData,substr(names(onePlateData),1,1))
    mcolname <- unique(substr(names(onePlateData),2,3))
    tmp <- lapply(tmp,function(x) {
        return( x[match(mcolname,substr(names(x),2,3))])
    })
    madata <- do.call(rbind,tmp)
    colnames(madata) <- mcolname
    rownames(madata) <- names(tmp)
    g <- pheatmap(madata,cluster_rows=FALSE,
             cluster_cols=FALSE,
             main=main,...)
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
#' @param genemap the mapping file of all samples
#' @param file the path of the file to generate to
#' @param outpath a character string naming the location the figures to
#' generate to
#' @param gDevice the graphic device (default: dev.new).
#' gDevice=NA generates no figure
#' @param ... arguments of the graphic device
#' @return a data frame with annotated information of each well
#' @import ggplot2
#' @export
#' @examples
#' data(demoCell)
#' data(genemap)
#' anno.data <- GenerateReport(list(oneCell),genemap,gDevice=NA)
#' colnames(anno.data)
#' 

GenerateReport <- function(lstCells,genemap,file=NULL,
                           outpath="./",gDevice="dev.new",...){
    if(any( !c("Barcode","Well") %in% colnames(genemap))){
        stop("The mapping table should contain Barcode and Well of each gene.")
    }
    genemap <- genemap[genemap$Well!="",]
    rownames(genemap) <- paste(genemap$Barcode,genemap$Well,sep=":")
    qc.data <- lstCells[[1]]@qc.data
    colnames(qc.data) <- paste(lstCells[[1]]["name"],
                               colnames(lstCells[[1]]["qc.data"]),sep=".")
    sig.data <- lstCells[[1]]["Sig"]$SigMat
    colnames(sig.data) <- paste(lstCells[[1]]["name"],
                                colnames(lstCells[[1]]["Sig"]$SigMat),sep=".")
    if(length(lstCells)>1){
        for(i in 2: length(lstCells)){
            qc.data.2 <-
              lstCells[[i]]["qc.data"][match(rownames(qc.data),
                                          rownames(lstCells[[i]]["qc.data"],))
                                    ,]
            colnames(qc.data.2) <- paste0(lstCells[[i]]["name"],
                                          colnames(lstCells[[i]]["qc.data"]))
            qc.data <- data.frame(qc.data,qc.data.2)
            sig.data.2 <-
              lstCells[[i]]["Sig"]$SigMat[match(rownames(sig.data),
                                             rownames(lstCells[[i]]["Sig"]$SigMat,)
                                             ),]
            colnames(sig.data.2) <- paste(lstCells[[i]]["name"],
                                          colnames(lstCells[[i]]["Sig"]$SigMat),
                                          sep=".")
            sig.data <- data.frame(sig.data,sig.data.2)
        }
    }
    qc.data <- qc.data[match(rownames(genemap),rownames(qc.data)),]
    sig.data <- as.matrix(sig.data[match(rownames(genemap),rownames(sig.data)),]
                          )
    mode(sig.data) <- "numeric"
    l3.data <- cbind(genemap,qc.data,sig.data)
    threshold <- unlist(lapply(lstCells,function(obj){
        th <- obj["Sig"]$threshold[c("Low","High")]
        return(th)
    }))
    stats <- unlist(lapply(lstCells,function(obj){
        st <- obj["Sig"]$stats[c("Low","High")]
        return(st)
    }))
    statMat <- cbind(c("threshold","number"),rbind(threshold,stats))
    if(!is.null(file)){
        write("Analysis Report",file=file)
        suppressWarnings(write.table(statMat, row.names=FALSE, sep="\t",
                                     quote=FALSE, file=file, append=TRUE))
        cat("\n",file=file,append=TRUE)
        suppressWarnings( write.table(l3.data, file=file, append=TRUE,
                                      sep="\t", quote=FALSE,
                                      row.names=FALSE, col.names=TRUE))

    }
    if(!is.na(gDevice)){
        if(gDevice=="dev.new"){
            dev.new(); 
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice)
                 )(file.path(outpath,paste("Statistics",gDevice,sep=".")),
                   ...)
        }
        count <- NULL; subType <- NULL; Type <- NULL
        ## Just for avoiding warnings in Rcheck
        stats.ggplot <- data.frame(count = stats,
                                   subType = names(stats),
                                   Type = rep(sapply(lstCells,
                                     function(obj)obj["name"]),each=2)
                                   )
        xlab.ggplot <- apply(stats.ggplot,1,function(x)
                             sub(paste0(x["Type"],"."),"",x["subType"]))
        print(ggplot(data=stats.ggplot, aes(x=subType, y=count)) +
              geom_bar(aes(fill=Type), stat="identity") +
              scale_x_discrete(labels = xlab.ggplot ) +
              geom_text(aes(label=count)) +
              xlab(""))
        if(gDevice !="dev.new"){
            graphics.off();
        }
    }
    return(l3.data)
}
