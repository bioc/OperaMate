### R code from vignette source 'OperaMate-vignette.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
if (requireNamespace("BiocStyle", quietly = TRUE)) {
    BiocStyle::latex()
}


###################################################
### code chunk number 2: demo_config
###################################################
file <- file.path(system.file("Test", package = "OperaMate"), "demoData", "demoParam.txt")
## file.show(file)


###################################################
### code chunk number 3: demo_genemap
###################################################
file <- file.path(system.file("Test", package = "OperaMate"), "demoData", "genemap.csv")
genemap <- read.csv(file)
head(genemap, n = 5)


###################################################
### code chunk number 4: run_pipeline
###################################################
library(OperaMate) 
configFile <- file.path(system.file("Test", package = "OperaMate"), "demoData", "demoParam.txt")
operaReport <- operaMate(configFile, gDevice = "png")
names(operaReport)


###################################################
### code chunk number 5: out_dir (eval = FALSE)
###################################################
## tempdir()


###################################################
### code chunk number 6: create_one_plate
###################################################
onePlate <- expData(name = "DSIMGA02-s1",
                     path = file.path( 
                       system.file("Test",package = "OperaMate"),
                       "Matrix", "130504-s1-02.txt" ),
                     rep.id = "s1", exp.id = "DSIMGA02", 
                     format = "Matrix")


###################################################
### code chunk number 7: create_one_cell
###################################################
oneCell <- cellData(name = "Average Intensity of Nuclei",
                    positive.ctr = c("H02", "J02", "L02"),
                    negative.ctr = c("C23", "E23", "G23"))
oneCell


###################################################
### code chunk number 8: import_data
###################################################
datapath <- file.path(system.file("Test", package = "OperaMate"), "Matrix")
lstPlates <- loadAll(cellformat = "Matrix", datapath = datapath,
                     egFilename <- list(eg.filename = "Tab.130504-s1-01.txt",
                                        rep.id = "s1", exp.id = "01", sep = "-",
                                        barcode = "DSIMGA01")
                     )
oneCell <- cellData(name = "Average Intensity of Nuclei")
oneCell <- cellLoad(oneCell, lstPlates, neglect.well = c("*02", "*23"),
                    positive.ctr = c("H02", "J02", "L02"),
                    negative.ctr = c("C23", "E23", "G23"))
str(oneCell["origin.data"])


###################################################
### code chunk number 9: norm_cell
###################################################
oneCell <- cellNorm(oneCell, norm.method = "MP") 
str(oneCell["norm.data"])


###################################################
### code chunk number 10: batch_methods
###################################################
cellViz(oneCell, data.type = c("raw", "norm"), plateID = 1:6, outpath = tempdir())


###################################################
### code chunk number 11: visualize_onePlate
###################################################
cellViz(oneCell, data.type = c("raw", "norm"), plateID = 1, outpath = tempdir())


###################################################
### code chunk number 12: load_cellnum
###################################################
cell.cellNum <- cellData(name = "Cells.Analyzed")
cell.cellNum <- cellLoad(cell.cellNum, lstPlates, neglect.well = c("*02", "*23"),
                    positive.ctr = c("H02", "J02", "L02"),
                    negative.ctr = c("C23", "E23", "G23"))
oneCell <- cellNumLoad(oneCell, cell.cellNum)


###################################################
### code chunk number 13: qc_cell
###################################################
oneCell <- cellQC(oneCell, qcType = c("plateCorrelation", "wellSd", "cellNumber"),
                  qc.threshold = c(correlation = 0.7), outpath = tempdir())
str(oneCell["qc.data"])
head(oneCell["plate.quality"])


###################################################
### code chunk number 14: sig_detection
###################################################
oneCell <- cellSig(oneCell, method = "stable", th = c(0.05, 0.05),
                   outpath = tempdir())
names(oneCell["Sig"])


###################################################
### code chunk number 15: sig_vcplot
###################################################
labels <- c("Axin1") 
names(labels) <- c("DSIMGA04:C07")
cellSigPlot(oneCell, highlight.label = labels, outpath = tempdir())


###################################################
### code chunk number 16: sig_analysis
###################################################
genemap <- read.csv(file.path(system.file("Test", package = "OperaMate"),
                              "demoData", "genemap.csv"), stringsAsFactors = FALSE)
chart <- cellSigAnalysis(oneCell, genemap, organism = "mmusculus")
head(chart, n = 5)
cellSigAnalysisPlot(chart, prefix = oneCell@name, outpath = tempdir())


###################################################
### code chunk number 17: report
###################################################
report <- generateReport(list(oneCell), genemap, verbose = FALSE,
                         plot = FALSE)
head(report, n = 5)


###################################################
### code chunk number 18: session_info
###################################################
sessionInfo() 


