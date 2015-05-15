#' Data extraction from one report
#'
#' Extracts data in the report to the slot \code{data} in the \code{expData}
#' object. An inner function of \code{\link{loadAll}}.
#'
#' @param onePlate an expData object
#' 
#' @return an \code{expData} object with initialized slot \code{data}.
#' @docType methods
#' @examples
#' data(platemap)
#' platemap$path <- file.path(system.file("Test",package = "OperaMate"),
#' platemap$path)
#' lstPlates <- loadAll(cellformat="Tab",datapath=dirname(platemap$path[1]))
#' onePlate <- lstPlates[[1]]
#' parseTemplete(onePlate)
#' @export
parseTemplete <- function(onePlate){
    if(onePlate["format"]=="Tab"){
        x <- parseTab(onePlate)
    }else if (onePlate["format"]=="Matrix"){
        x <- parseMatrix(onePlate)
    }else{
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
parseTab <- function(onePlate){
    if(onePlate["format"]!="Tab"){
        stop("Error processor for ",onePlate["name"],": ",onePlate["format"],".")
    }
    expdata <- read.delim(onePlate["path"],sep="\t",stringsAsFactors=FALSE,skip=3)
    if(colnames(expdata)[1]!="Row"){
        stop("Error processor for ",onePlate["name"],": ",onePlate["format"],".")
    }
    WellName <- paste0(expdata[,"Row"],
                       formatC(expdata[,"Column"],width=2,flag="0"))
    twellID <- ( !is.na(expdata[,3]) & expdata[,3]!="" )
    expdata <- expdata[ twellID,]
    dataLoad( onePlate, as.list(expdata[,4:ncol(expdata)]), WellName[twellID] )
}


## Data extraction when the report is of matrix format
##
## Extracts data in the report which is of matrix format. An inner function
## of \code{parseTemplete}.
##
## @param onePlate an expData object
## @return an \code{expData} object with initialized slot \code{data}.
## @docType methods
parseMatrix <- function(onePlate){
    if(onePlate["format"]!="Matrix"){
        stop("Error processor for ",onePlate["name"],": ",onePlate["format"],".")
    }
    expdata <- read.delim(onePlate["path"],
                          sep="\t",stringsAsFactors=FALSE,skip=3,header=FALSE)
    sline <- grep("Parameter: ",expdata[,1])
    ## Format Checking
    region <- expdata[(sline[1]+2):(sline[2]-2),2:ncol(expdata)]
    rownames(region) <- expdata[(sline[1]+2):(sline[2]-2),1]
    colnames(region) <- formatC(seq(1,ncol(region)),width=2,flag="0")
    vecregion <- as.vector(as.matrix(region))
    WellName <- apply(expand.grid(rownames(region),colnames(region))
                      ,1,paste,collapse="")
    twellID <- 
      (!is.na(as.vector(as.matrix(region)))&as.vector(as.matrix(region))!="")
    lstmatrix <- lapply(sline[2:length(sline)],function(line){
        fmatrix <-
          expdata[(line+2):(line+1+nrow(region)),2:ncol(expdata)]
        f.3m <- as.vector(as.matrix(fmatrix))
        mode(f.3m) <- "numeric"
        return(f.3m[twellID])
    })
    WellName <- WellName[twellID]
    names(lstmatrix) <- gsub(" ",".",expdata[sline[2:length(sline)],2])
    dataLoad( onePlate, lstmatrix, WellName)
}   

#' Data importing
#' 
#' Initializes a list of \code{expData} objects from the Columbus system
#' reports.
#'
#' "Tab" and "Matrix" are corresponding to table and matrix format in the
#' Columbus system report. They requires all reports of the same format, in
#' the same location, and with the file names as  *-replicateID-plateID.txt.
#' Otherwise,choose "Others" and specify the parameter of each file in
#' \code{platemap}. If the reports are of other cellformats, you can specify
#' its cellformat and rewrite the function \code{parseTemplete} to import the
#' data seperately. See the example to find the format for the slot \code{data}
#' and \code{wellID}.
#' 
#' @param cellformat character specifying the format of the reports
#' @param datapath character specifying the location of the reports
#' @param platemap data frame, required when cellformat = "Others".
#' See an example as \code{\link{platemap}}.
#' @param prefix character, the Barcode prefix, e.g. "DSIMGA"
#' 
#' @return a list of \code{expData} objects
#' @docType methods
#' @examples
#' data(platemap)
#' platemap$path <- file.path(system.file("Test",package = "OperaMate"),
#' platemap$path)
#' ## Not consistent format
#' lstPlates <- loadAll(cellformat="Others",platemap=platemap)
#' ## Consistent format
#' ## lstPlates <- loadAll(cellformat="Tab",datapath=dirname(platemap$path[1]))
#' ## User defined format
#' ## lstPlates <- loadAll(dirname(platemap$path[1]),cellformat="usrDef")
#' ## lstPlates <- lapply(lstPlates,function(onePlate) parseTemplete(onePlate))
#' lstPlates[[1]]
#' @export
loadAll <- function(cellformat,datapath="./",platemap=NULL, prefix = "DSIMGA" ){
    if ( cellformat == "Others" ){
        if( is.null(platemap) ){
            stop("Please specify the generating data formats for all 
                 experiments if they are not unique.
                 Try data(platemap) for the required format of the
specification file.")
        }
        ## Check if the existence and legality of basic information
        if( any( !c("file.name","format.type","Barcode","rep.id")
                %in% colnames(platemap) ) ){
            stop("The format of mapping file is illegal!
Try data(platemap.example) for the required format of the specification file.")
        }
        if( any( !platemap$format.type %in% c("Tab","Matrix" )) ){
            warning("Please choose a proper format from ( Tab, Matrix ).
Otherwise, you need to redefine the function parseTemplete to import data.")
        }
        lstPlates <- apply(platemap,1,function(x){
            onePlate <-
              expData( name  = paste(x["Barcode"],x["rep.id"],sep="-"),
                      path   = x["path"],
                      rep.id = x["rep.id"],
                      exp.id = x["Barcode"],
                      format = x["format.type"])
            return( parseTemplete(onePlate) )
        })
    }else{
        ## File names are in a standard way
        ## Ex. 20140101-s1-02.txt (date-replication-plateID)
        if (!cellformat %in% c("Tab","Matrix")){
            warning("Please choose a proper format from ( Tab, Matrix ).
Otherwise, you need to redefine the function parseTemplete to import data.")
        }
        pfiles <- list.files( datapath,recursive=TRUE )
        lstPlates <- lapply(pfiles,function(x){
            pinfo <- unlist(strsplit(
              sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
              ,"-"))
            exp.id <- paste0(prefix,pinfo[3])
            onePlate <-
              expData(  name = paste(exp.id,pinfo[2],sep="-"),
                        path = file.path(datapath,x),
                      exp.id = exp.id,
                      rep.id = pinfo[2],
                      format = cellformat)
            return( parseTemplete(onePlate) )
        })
    }
    return(lstPlates)
}

