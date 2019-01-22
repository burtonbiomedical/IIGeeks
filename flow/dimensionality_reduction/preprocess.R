#### LOAD IN FLOW DATA, TRANSFORM, AND COMPENSATE ####

# Dependencies
require(openCyto)
require(ggcyto)
require(MASS)

#Load in data
fcsFiles <- dir("/home/rossco/Documents/IIGeeks/flow/fcs", 
                full.names = TRUE, include.dirs=FALSE)
fcsFiles <- fcsFiles[!file.info(fcsFiles)$isdir]
ncfs  <- read.ncdfFlowSet(fcsFiles)
gs <- GatingSet(ncfs)

#Compensate data
comp_matrix <- ncfs[[1]]@description$SPILL
gs <- compensate(gs, comp_matrix)

#Transform data
trans <- estimateLogicle(gs[[1]], channels = colnames(comp_matrix))
gs <- transform(gs, trans)

#Remove Debris
chnls <- c("FSC-A", "SSC-A")
nonDebris <- rectangleGate(filterId = "nonDebris", list("FSC-A" = c(35000, 125000),
                                                        "SSC-A" = c(10000,50000)))
add(gs, nonDebris, parent="root")
recompute(gs)

#Get Singlets
fr <- getData(gs, "nonDebris")[[1]]
singlet_gate <- openCyto:::.singletGate(fr, channels=c("FSC-A","FSC-H"), 
                                        wider_gate=TRUE, prediction_level=0.99, maxit=15)
add(gs, singlet_gate, parent="nonDebris", name="single cells")
recompute(gs)

#Remove dead cells
fr <- getData(gs, "single cells")[[1]]
liveCD3 <- openCyto::flowClust.2d(fr, "PE-Cy5-A","AmCyan-A", K=2, quantile = 0.9)
add(gs, liveCD3, parent="single cells", name="live CD3")
recompute(gs)

#Get the intensity data
fr <- getData(gs, "liveCD3")[[1]]
m <- exprs(fr)

#Give meaningful column names
meaningfulNames <- function(fr, m){
  c_names <- colnames(fr)
  desc <- fr@parameters$desc
  for(i in seq(1,length(c_names))){
    if(!is.na(desc[i])){
      c_names[i] <- desc[i]
    }
  }
  colnames(m) <- paste(c_names)
}

meaningfulNames(fr, m)

#Write to disc
write.matrix(m, file="FULL_PANEL.csv", sep=",")
