#' @import stats 
#' @import methods
#' @importFrom MASS fitdistr
NULL

#' Batch effect visualization
#' 
#' Visualizes batch effects by three ways: hieratical clustering (heatmap),
#' the principal component analysis (PCA) and the median scatter plot (median).
#'
#' Heatmap: a large region of distinguishing color indicates the systematical
#' biases; PCA: points distant from others indicate data with excessive noises;
#' Median: the distribution of the plate and well medians.
#'
#' @param object a cellData object
#' @param exps exps = c("all", "exp","ctr"), "exp" and "ctr" consider only the
#' sample and control wells respectively. 
#' @param type type = c("heatmap","median","PCA")
#' @param ifCompare logic, comparison before and after normalization will be
#' shown if the data normalization has been processed.
#' @param outpath a character string naming the location the figures to generate
#' to
#' @param gDevice the graphic device (default: dev.new)
#' @param control.par list, the parameter to be passed to plotting methods
#' @param ... arguments of the graphic device
#' @section control.par:
#' type="heatmap"
#' \describe{
#' \item{\code{    fontsize}:}{the fontsize of the plot (default: 15)}
#' \item{\code{    color}:}{vector of character specifying the colors used in
#' heatmap}
#' \item{\code{    fontsize_row}:}{fontsize of rownames (default: 0)}
#' \item{\code{    fontsize_col}:}{fontsize of colnames (default: 0)}
#' \item{\code{    xlab}:}{character, X axis label}
#' \item{\code{    ylab}:}{character, Y axis label}
#' \item{\code{    remove.noise}:}{logical. If TRUE, the heatmap only represents
#' data in the range of [median-3mad, median+3mad]}
#' \item{\code{    ...}:}{arguments to be passed to \code{pheatmap}}
#' }
#'
#' type="median"
#' \describe{
#' \item{\code{    lwd}:}{the line widths}
#' \item{\code{    pch}:}{the plotting symbols}
#' \item{\code{    main}:}{the title of the figure}
#' \item{\code{    ...}:}{other parameters for function \code{plot}}
#'}
#' 
#' type="PCA"
#' \describe{
#' \item{\code{    main}:}{the title of the figure}
#' \item{\code{    pch}:}{the plotting symbols of the legends}
#' \item{\code{    ...}:}{other parameters for function \code{plot}}
#' }
#' 
#' @examples
#' 
#' data(demoCell)
#' cellViz(oneCell,exps="exp",type="heatmap",ifCompare=FALSE)
#' cellViz(oneCell,exps="exp",type="median",ifCompare=FALSE)
#' cellViz(oneCell,exps="exp",type="PCA")
#' @return Invisibly a list of two components:
#' \describe{
#' \item{\code{g1}:}{Plotting results before data normalization.}
#' \item{\code{g2}:}{Plotting results after data normalization.
#' g2=NULL if ifCompare=FALSE.}
#' }
#' 
#' @docType methods
#' @export
#' 
cellViz <- function(object,exps=c("all","exp","ctr"),
                    type=c("heatmap","median","PCA"),gDevice="dev.new",
                    outpath="./", ifCompare=TRUE,control.par=list(),...){
    method.par <- list()
    g1 <- NULL; g2 <- NULL
    if(type=="heatmap"){
        method.par <- list(fontsize=15,
                           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                           fontsize_row=0, fontsize_col=0,xlab="WellID",ylab="PlateID",
                           remove.noise=TRUE)
    }else if (type == "median"){
        method.par <- list(main="",cex.main=1,lwd=1.5,pch=16, type="b")
    }else if (type == "PCA"){
        method.par <- list(main="",cex.main=1,lwd=1.5,pch=19,
                           color.background="gray")
    }
    control.par1 <- control.par[names(control.par)%in%names(method.par)]
    method.par[names(control.par1)] <- control.par1

    if(!length(object["norm.data"])){
        ifCompare <- FALSE
    }
    showID <- rownames(object["origin.data"])
    if (exps=="exp"){
        showID <- object["expwell"]
    }else if (exps=="ctr"){
        showID <- object["ctrwell"]
    }
    if( type == "heatmap" ){
        if(gDevice=="dev.new"){
            dev.new();
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice))(
                              file.path(outpath,paste0(type,".",object["name"],
                                                       ".before.normalization.",gDevice)),...)
        }
        g1 <- cellHeatmap(object["origin.data"][showID,],fontsize=method.par$fontsize,
                    color=method.par$color,
                    fontsize_row=method.par$fontsize_row,
                    fontsize_col=method.par$fontsize_col,
                    xlab=method.par$xlab,
                    ylab=method.par$ylab,
                    remove.noise=method.par$remove.noise)
        if(gDevice !="dev.new"){
            graphics.off();
        }
        if(ifCompare){
            if(gDevice=="dev.new"){
                dev.new();
            }else{
                graphics.off()
                eval(expr=parse(text=gDevice))(
                                  file.path(outpath,paste0(type,".",object["name"],
                                                           ".after.normalization.",gDevice)),...)
            }
            g2 <- cellHeatmap(object["norm.data"][showID,],fontsize=method.par$fontsize,
                        color=method.par$color,
                        fontsize_row=method.par$fontsize_row,
                        fontsize_col=method.par$fontsize_col,
                        xlab=method.par$xlab,
                        ylab=method.par$ylab,
                        remove.noise=method.par$remove.noise)
            if(gDevice !="dev.new"){
                graphics.off();
            }
        }
    }else{
        if(gDevice=="dev.new"){
            dev.new(); 
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice))(
                              file.path(outpath,paste0(type,".",object["name"],".",gDevice)),...)
        }
        if ( type == "median" ){
            if(ifCompare){
                par(mfrow=c(2,2))
            }else{
                par(mfrow=c(1,2))
            }
            g1 <- cellMedian(object["origin.data"][showID,],main="Before Normalization",
                       lwd=method.par$lwd,pch=method.par$pch)
        }else if(type=="PCA"){
            if(ifCompare){
                par(mfrow=c(2,2))
            }else{
                par(mfrow=c(1,2))
            }
            g1 <- cellPCA(object["origin.data"][showID,],
                    main="Before Normalization",pch=method.par$pch)
        }
        if(ifCompare){
            if ( type == "median" ){
                g2 <- cellMedian(object["norm.data"][showID,],main="After Normalization",
                           lwd=method.par$lwd,pch=method.par$pch)
            }else if(type=="PCA"){
                g2 <- cellPCA(object["norm.data"][showID,],main="After Normalization",
                        pch=method.par$pch)
            }
        }
        if(gDevice !="dev.new"){
            graphics.off();
        }
    }
    invisible(list(g1=g1,g2=g2))
}



## Heatmap of data matrix
## 
## Performs hieratical clustering to data matrix.
## An inner function of \code{\link{cellViz}}
## 
## @param data the data matrix
## @param fontsize the fontsize of the plot (default: 15)
## @param color vector of character specifying the colors used in heatmap
## @param fontsize_row fontsize of rownames (default: 0)
## @param fontsize_col fontsize of colnames (default: 0)
## @param xlab character, X axis label
## @param ylab character, Y axis label
## @param remove.noise logical. If TRUE, the heatmap only represents data in
## the range of [median-3mad, median+3mad]
## @param ... arguments to be passed to \code{pheatmap}
## @return Invisibly a list of component by \code{pheatmap}.
#' @importFrom grid pushViewport grid.newpage grid.text gpar
## @docType methods
#' @import pheatmap
cellHeatmap <- function(
  data, fontsize=15,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row=0, fontsize_col=0,xlab="WellID",
  ylab="PlateID",remove.noise=TRUE,...){
    setHook( "grid.newpage", function()
            pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp",
                                  just=c("right","top"))), action="replace" )
    if(remove.noise){
        medianData <- median(as.vector(data),na.rm=TRUE)
        madData <- mad(as.vector(data),na.rm=TRUE)
        breaks <- seq(medianData-3*madData,
                      medianData+3*madData,length.out=51)
    }else{
        breaks <- NA
    }
    g <- pheatmap( data,color=color,
                  fontsize=fontsize,fontsize_row=fontsize_row,
                  fontsize_col=fontsize_col,
                  breaks = breaks,...)
    setHook( "grid.newpage", NULL, "replace" )
    grid.text( xlab, x=-0.02, rot=90, gp=gpar(fontsize=fontsize) )
    grid.text( ylab, y=-0.01, gp=gpar(fontsize=fontsize) )
    invisible(g)
}

## The distributions of plate and well medians
##
## Plots the distribution of plates and well medians.
## An inner function of \code{\link{cellViz}}.
## @param data the data matrix
## @param lwd the line widths
## @param pch the plotting symbols
## @param main the title of the figure
## @param ... other parameters for the function \code{plot}
## @return Invisibly a list with two components:
## \describe{
## \item{\code{medCol}:}{the median of all wells in each plate}
## \item{\code{medRow}:}{the median of each well location across all plates}
## }
## @docType methods
## 
cellMedian <- function(data,lwd=1.5,pch=16,main="",...){
    lstmedianCol <- apply(data,2,median,na.rm=TRUE)
    plot( 1:length(lstmedianCol), sort(lstmedianCol),
         type="b", lwd=lwd, pch=pch,
         ylab="Median intensity of Plate",
         xlab="PlateID",
         main=paste("Plate Median Distribution",main),...)
    lstmedianRow <- apply(data,1,median,na.rm=TRUE) 
    plot( 1:length(lstmedianRow), sort(lstmedianRow),
         type="b", lwd=lwd, pch=pch,
         ylab="Median intensity of Well",
         xlab="WellID",
         main=paste("Well Median Distribution",main),...)
    invisible(list(medCol=lstmedianCol, medRow=lstmedianRow))
}

## Batch effect visualization by principal component analysis (PCA)
##
## Represents data by their two pricipal components.
## An inner function of \code{\link{cellViz}}.
##
## Comparison before and after normalization will be shown
## if the data normalization has been processed.
## 
## @param data the data matrix
## @param pch the plotting symbols of the legends
## @param main the title of the figure
## @param ... arguments to be passed to \code{plot}
## @return Invisibly a list with class "princomp" by \code{princomp}
## @docType methods
## 
cellPCA <- function(data, main="", pch=19,...){
    data[is.na(data)] <- mean(data,na.rm=TRUE)
    if((ncol(data)/2)>nrow(data)){
        warning("'princomp' can only be used with more units than variables")
        return()
    }else if( ncol(data)>nrow(data) ){
        spN <- sample(colnames(data),nrow(data))
    }else{
        spN <- colnames(data)
    }
    PcData <- princomp(data[,spN],cor=TRUE)
    if(nchar(rownames(data)[1])==3){
        tag1 <- substr(rownames(data),1,1)
        tag2 <- substr(rownames(data),2,3)
    }else if (nchar(rownames(data)[1])==4){
        tag1 <- substr(rownames(data),1,2)
        tag2 <- substr(rownames(data),3,4)
    }
    lstColor1 <- split(rownames(data),tag1)
    lgTag1 <- names(lstColor1)
    names(lstColor1) <- rainbow(n=length(lstColor1))
    color1 <- unlist(lapply(names(lstColor1),function(x){
        ctmp <- rep(x,length(lstColor1[[x]]))
        names(ctmp) <- lstColor1[[x]]
        return(ctmp)
    }))
    color1 <- color1[match(rownames(data),names(color1))]
    lstColor2 <- split(rownames(data),tag2)
    lgTag2 <- names(lstColor2)
    names(lstColor2) <- rainbow(n=length(lstColor2))
    color2 <- unlist(lapply(names(lstColor2),function(x){
        ctmp <- rep(x,length(lstColor2[[x]]))
        names(ctmp) <- lstColor2[[x]]
        return(ctmp)
    }))
    color2 <- color2[match(rownames(data),names(color2))]
    plot(PcData$scores[,1],PcData$scores[,2],
         col=color1,type="p",
         xlab = paste("PC1(",PcData$sdev[1],")",sep=""),
         ylab = paste("PC2(",PcData$sdev[2],")",sep=""),
         main = main,
         pch=pch,...)
    legend("bottomright",legend=lgTag1,col=names(lstColor1), pch=pch)
    plot(PcData$scores[,1],PcData$scores[,2],
         col=color2,type="p",
         xlab = paste("PC1(",PcData$sdev[1],")",sep=""),
         ylab = paste("PC2(",PcData$sdev[2],")",sep=""),
         main = main,
         pch=pch,...)
    legend("bottomright",legend=lgTag2,col=names(lstColor2), pch=pch)
    invisible(PcData)
}


#' One plate visualization
#' 
#' Visualizes data of one plate using heatmap.
#' 
#' @param object a cellData object
#' @param ID character, the plate ID to be visualized
#' @param outpath a character string naming the location the figure to
#' generate to
#' @param gDevice the graphic device (default: dev.new),
#' gDevice=NA generates no figure
#' @param ... arguments to be passed to \code{pheatmap}
#' @docType methods
#' @import pheatmap
#' @export
#' @return Invisibly a list of two components:
#' \describe{
#' \item{\code{g1}:}{a list of component by \code{pheatmap} based on data
#' before normalization}
#' \item{\code{g2}:}{a list of component by \code{pheatmap} based on data
#' after normalization}
#' }
#' @examples
#' data(demoCell)
#' plateViz(oneCell,ID="DSIMGA04-s2")
#' 

plateViz <- function(object,ID,gDevice="dev.new",outpath="./",...){
    data <- object["origin.data"][,ID]
    g1 <- NULL; g2 <- NULL
    if(gDevice=="dev.new"){
        dev.new();
    }else{
        graphics.off()
        eval(expr=parse(text=gDevice)
             )( file.path(outpath,paste0(ID, ".before.normalization.",gDevice))
               )
    }
    g1 <- showPlate(data,main=paste(ID,"Before Normalization"),...)
    if(gDevice!="dev.new"){
        graphics.off();
    }
    if(length(object["norm.data"])){
        data <- object["norm.data"][,ID]
        if(gDevice=="dev.new"){
            dev.new();
        }else{
            graphics.off()
            eval(expr=parse(text=gDevice)
                 )(file.path(outpath,paste0(ID,".after.normalization.",gDevice))
                   )
        }
        g2 <- showPlate(data,main=paste(ID,"After Normalization"),...)
        if(gDevice!="dev.new"){
            graphics.off();
        }
    }
    invisible(list(g1=g1,g2=g2))
}



#' Hits volcano plot
#' 
#' Visualizes hits by volcano plot.
#'
#' Users can highlight a certain samples during plottting.
#'
#' @param object a cellData object
#' @param outpath a character string naming the location the figure to
#' generate to
#' @param gDevice the graphic device (default: dev.new)
#' @param color.highlight a character specifying the color of the hits
#' @param color.background a character specifying the color of the other samples
#' @param highlight.label a vector of characters specifying the names of the
#' samples to be highlighted, with the names are the "barcode:wellID".
#' @param ... arguments of the graphic device
#' @docType methods
#' @examples
#' data(demoCell)
#' labels <- c("Axin1")
#' names(labels) <- c("DSIMGA04:C07")
#' cellSigplot(oneCell,highligh.label=labels)
#' @import ggplot2
#' @export
#' @return  Invisibly an object of \code{ggplot}
#' 
cellSigplot <- function(object,outpath="./",gDevice="dev.new",
                        color.highlight = "red",
                        color.background = "blue",
                        highlight.label = NULL,
                        ...
                        ){
    
    if(is.null(object["Sig"]$pvalue)){
        stop("This method requires three or more replicated data.")
    }
    objSig <- object["Sig"]
    objData <- object["qc.data"]
    vc.pvalue <- -log10(objSig$pvalue)
    vc.Data <- log2(objData[,"mean"])
    vc.threshold <- ifelse(apply(objSig$SigMat,1,any),
                           color.highlight,color.background)
    vc.threshold[is.na(vc.threshold)] <- color.background
    volcanoData <- data.frame(exp = vc.Data ,
                              pvalue = vc.pvalue,
                              threshold = vc.threshold)
    if(gDevice=="dev.new"){
        dev.new()
    }else{
        graphics.off()
        eval(expr=parse(text=gDevice)
             )(file.path(outpath,paste0(object["name"],".vocanoPlot.",gDevice)),...
               )
    }
    
    exp <- NULL; pvalue <- NULL; label <- NULL
    ## Just for avoiding warnings in Rcheck
    g = ggplot( data=volcanoData, aes(x=exp, y=pvalue))+
      geom_point( size=2.5,shape=21,colour="black",
                 fill=volcanoData$threshold,alpha=1) +
                   theme(legend.position = "none") +
                     xlab("log2( intensity )") +
                       ylab("-log10 p-value")
    if( !is.null(highlight.label) ){
        TextFrame <- volcanoData[match(names(highlight.label),
                                       rownames(volcanoData)),]
        if(is.null(TextFrame)){
            warning("Cannot find the highlight well IDs.
Parameter 'highlight.label' is a vector of labels with the vector names are in
the format of plateID:wellID.")
        }else{
            TextFrame$label <- highlight.label
            TextFrame <- transform(TextFrame,
                                   w=sapply(label,strwidth,"inches"),
                                   h=sapply(label,strheight,"inches"))
            w <- NULL;  h <- NULL      ## Just for avoiding warnings in Rcheck
            g <- g +
              geom_point(alpha=1,data=TextFrame,aes(x=exp,y=pvalue),
                         fill="yellow", colour="black",shape=21,size=3)+
                           geom_rect(data = TextFrame,
                                     aes(xmin = exp -w/2,
                                         xmax = exp + w/2,
                                         ymin = pvalue + h/2,
                                         ymax = pvalue + 3*h/2),
                                     color="yellow",fill = "yellow") +
                                       geom_text(data=TextFrame,
                                                 aes(x=exp,y=pvalue,
                                                     label = label),
                                                 colour="black",size=4,vjust=-1)
        }
    }
    suppressWarnings(print (g))
    if(gDevice!="dev.new"){
        graphics.off();
    }
    invisible(g)
}


#' Hits function analysis
#'
#' Performs function analysis using DAVID Web Service
#'
#' @param object a cellData object
#' @param genemap the mapping file of all samples
#' @param email the email for the DAVID WebService registration
#' @param type type=c("all","high","low"), all hits or only the high/low
#' expressed hits to analysis
#' @param david.terms a list of DAVID annotation categories
#'
#'                    Default: list(GO=c("GOTERM_BP_FAT","GOTERM_MF_FAT",
#' "GOTERM_CC_FAT"),PATHWAY=c( "KEGG_PATHWAY","REACTOME_PATHWAY","BIOCARTA"))
#' @param count numeric, number of items
#' @param threshold numeric, threshold in the DAVID analysis
#' @param outpath a character string naming the location the figures to
#' generate to
#' @param gDevice the graphic device (default: dev.new).
#' gDevice=NA generates no figure
#' @param showTermNum numeric, number of terms to show in the figure
#' @param showTermBy character, the order of the DAVID chart,
#' showTermBy is one of the column names of the chart.
#' @param showTermCol character, color of the bars in the figure
#' @param file the path of the DAVID report to generate to.
#' file=NA disables report generation.
#' @param ... arguments of the graphic device
#' @return a data frame of the DAVID report with annotated well IDs
#' and gene symbols
#' @docType methods
#' @examples
#' data(demoCell)
#' data(genemap)
#' ## Please register your email at DAVID WebService
#' chart <- cellSigAnalysis(oneCell,genemap,
#'                         email="cliu@@sjtu.edu.cn",type="all",
#'                         file=NA)
#' colnames(chart)
#' 
#' @importFrom RDAVIDWebService DAVIDWebService addList setAnnotationCategories
#' @importFrom RDAVIDWebService getFunctionalAnnotationChart  is.connected
#' @export
#'

cellSigAnalysis <- function(object,genemap,email,
                            type=c("all","high","low"),
                            david.terms = NULL,count=2L,threshold=0.1,
                            showTermNum=NA,showTermBy="PValue",
                            showTermCol="blue",file=NA,outpath="./",
                            gDevice="dev.new",...){
    if(missing(email)|| is.na(email)){
        stop("Please register according to the instruction of DAVID Web Service
at http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.")
    }
    if(is.null(david.terms)){
        david.terms <- list(
          GO=c("GOTERM_BP_FAT","GOTERM_MF_FAT", "GOTERM_CC_FAT"),
          PATHWAY=c( "KEGG_PATHWAY","REACTOME_PATHWAY","BIOCARTA"))
    }
    colnames(genemap)[colnames(genemap)=="Gene.ID"] <- "GENE.ID"
    if(any( !c("Barcode","Well","GENE.ID") %in% colnames(genemap))){
        stop("The mapping table should contain barcode and Well of each gene and
 the Entrez Gene ID.")
    }
    objSig <- object["Sig"]
    sig <- objSig$SigMat[match(paste(genemap$Barcode,genemap$Well,sep=":"),
                                   rownames(objSig$SigMat)),]
    if(type=="high"){
        ID <- rownames(objSig$SigMat)[objSig$SigMat[,"High"]]
    }else if(type == "low"){
        ID <- rownames(objSig$SigMat)[objSig$SigMat[,"Low"]==1]
    }else{
        ID <- rownames(objSig$SigMat)[objSig$SigMat[,"Low"]==1
                                          |objSig$SigMat[,"High"]==1]
    }
    genes <- as.character(
      genemap[match(na.omit(ID),paste(genemap$Barcode,genemap$Well,sep=":")),
              "GENE.ID"])
    david <- try(DAVIDWebService$new(email=email),silent=TRUE)
    if(class(david)=="try-error"){
        warning("Cannot connect to DAVID Webservice")
        return(NULL)
    }
    addList(david,genes,idType="ENTREZ_GENE_ID",
            listName=paste(object["name"],type[1],sep="."),listType="Gene")
    setAnnotationCategories(david, unlist(david.terms))
    chart <- getFunctionalAnnotationChart(david,threshold=threshold,count=count)
    wellID <- sapply( strsplit(chart$Genes,","),function(x){
        paste(apply(genemap[match(as.numeric(x),genemap$GENE.ID),
                            c("Barcode","Well")]
                    ,1,paste,collapse="-"),collapse=";")})
    geneSymbol <- sapply( strsplit(chart$Genes,","),function(x){
        paste(genemap[match(as.numeric(x),genemap$GENE.ID),"GeneSymbol"]
              ,collapse=";")})
    chart.new <- data.frame(chart,wellID,geneSymbol)
    if(!is.na(file)){
        write.table(chart.new,sep="\t",quote=FALSE,row.names=FALSE,
                    file=file)
    }
    if(!is.na(gDevice)){
        for(term in names(david.terms)){
            chart0 <- chart[chart$Category %in% david.terms[[term]],]
            chart0 <- chart0[order(chart0[,showTermBy],decreasing=TRUE),]
            if(nrow(chart0)==0) {
                next
            }
            if(!is.na(showTermNum)){
                chart0 <- tail(chart0,n=showTermNum)
            }
            cellSigAnalysisPlot(chart0,prefix=paste(object["name"],type[1],sep=".")
                                ,type=term,gDevice=gDevice,outpath=outpath,
                                col=showTermCol,...)
        }
    }
    return(chart.new)
}
## DAVID annotation visualization
##
## Visualizes the DAVID annotation chart by barplot.
##
## @param chart data frame, the functional annotation chart
## @param outpath a character string naming the location the figures to
## generate to
## @param gDevice the graphic device (default: dev.new).
## gDevice=NA generates no figure
## @param prefix character name of the figure
## @param type character, the category of the chart to show
## @param col character, color of the bars in the figure
## @param ... arguments of the graphic device
## @return Invisibly a numeric vector by \code{barplot}
## 
## @docType methods

cellSigAnalysisPlot <- function(chart,prefix="Sig",type=NULL,
                                gDevice="dev.new", outpath="./",col="blue",...){
    if(is.null(type)){
        type <- chart[1,"Category"]
    }else if (length(type)>1){
        type <- paste(type,collapse=".")
    }
    if(gDevice=="dev.new"){
        dev.new();
    }else{
        graphics.off()
        eval(expr=parse(text=gDevice)
             )(file.path(outpath,paste0(prefix,"-",type,".",gDevice)),...)
    }
    pvalue <- (-log10(chart$PValue))
    g <- barplot(pvalue,xlim=c(0,6),col=col,horiz=TRUE, names.arg=chart$Term,
            legend.text = type, args.legend = list(x = "bottomright"),
            xlab="-log10 p-value",cex.name=0.7)
    if(gDevice!="dev.new"){
        graphics.off();
    }
    invisible(g)
}



