#' @import grDevices
.onLoad <- function(libname, pkgname) {
#    assign("global_opts", new.env(), envir = parent.env(environment()))
#    global.opts <- options()
    if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
    assign("op", options(), parent.env(environment()))  #save old options

    options(list( stringsAsFactors = FALSE
                 ))
    op.operaMate <- list(
      opm.outpath = "./",
      opm.normalization.method = c("MP", "PMed", "Z", "Ctr", "None"),
      opm.QC.type = c("plateCorrelation", "wellSd", "zFactor", "cellNumber"),
      opm.grDevices = c("jpeg","png", "tiff","bmp", "pdf"),
      opm.QC.threshold = c(correlation = 0.8, zfactor = 0.5, cellnumber = 50),
      opm.replace.badPlateData = TRUE,
      opm.heatmap.color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
      opm.heatmap.fontsize = 15,
      opm.heatmap.fontsize_row = 10,   
      opm.heatmap.fontsize_col = 10,   
      opm.heatmap.rmnoise = TRUE,
      opm.adjust.methods = "fdr",
      opm.threshold.digits = 3,
      opm.sig.color.highlight = "red",
      opm.sig.color.background = "blue",
      opm.sig.label.color = "yellow",
      opm.filename.example = list(eg.filename = "0205-s2-01.txt", rep.id = "s2",
        exp.id = "01", sep = "-", barcode = "DSIMGA01")
      )
    toset <-!(names(op.operaMate) %in% names(op))
    if(any(toset)) options(op.operaMate[toset])
    invisible()
}

.onUnload <- function(libname, pkgname) {
    on.exit(options(op), add = TRUE)
}
