par(mfrow = c(1,3), mar = c(0,0,1,0))

xlabel <- "Genes"
phi <- 35; theta <- 45
feature <- "CDS"
source("3dplots-plot.r")

phi <- 48; theta <- 120
feature <- "INT"
source("3dplots-plot.r")

xlabel <- "Elements"
phi <- 35; theta <- 120
feature <- "TE"
source("3dplots-plot.r")
