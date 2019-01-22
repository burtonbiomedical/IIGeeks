require(FlowSOM)

fcsFiles <- dir("/home/rossco/Documents/IIGeeks/flow/fcs", 
                full.names = TRUE, include.dirs=FALSE)
fcsFiles <- fcsFiles[!file.info(fcsFiles)$isdir]

fSOM <- FlowSOM(fcsFiles, compensate = TRUE,transform = TRUE,toTransform=c(6:14),
                scale = TRUE, colsToUse = c(6:14), xdim = 7, ydim = 7,
                seed = 42, maxMeta = 15)

PlotStars(fSOM$FlowSOM, backgroundValues = as.factor(fSOM$metaclustering), view="grid")
