library(circlize)
ds <- readRDS("enrich-methyl-CLUSTERS.data")
links <- readRDS("genes_common_clusters.data")
GOs <- c("BP","MF","CC")
################################################################################
# INITIALISATION OF THE PLOT
# Initmat contains each compartment, repeated twice, with the limits to be
# used for each sector horizontal limit:
initmat <- data.frame(cluster = rep(unique(ds$cluster), each = 2),
                      width = c(0,3))

# Factor to maintain order:
fa <- unique(ds$cluster)
fa <- factor(initmat$cluster, levels = fa)
# Plot initialisation:
par(mar = c(0.1, 2.5, 0.1, 0.1))
circos.par("track.height" = 0.14)
circos.initialize(factors = fa, x = initmat$width)
################################################################################
# Outside region with names of development stages:
circos.trackPlotRegion(factors = fa, ylim = c(0, 1), track.height = 0.1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                       circos.text(1.5, 0.2,
                           paste("Cluster",get.cell.meta.data("sector.index")),
                                   niceFacing = TRUE) } )
################################################################################
# New inside region with histograms:
circos.trackPlotRegion( factors = fa, ylim = c(0,max(ds)) )

# HISTOGRAMS AND GO TEXT:
# Quick function to rescale each value vector to 0-1 range:
"reScale" <- function(x, from, to) (x-min(x)) / max(x-min(x)) * (to-from) + from

for (i in 1:max(initmat$cluster))           # For each cluster:
{
 # Move to that cluster:
 circos.updatePlotRegion(sector.index = i)
 ylim <- get.cell.meta.data("ylim")
 # Extract the scores that will be the rectangles' height:
 rectval <- ds[ds$cluster == i , 2:4]
 # Names and colours of the three bars:
 textval <- c("BP","MF","CC")
 colval <- c("#0066CC","#CCFFFF","#99CC99")

 for (j in 1:3) {
    circos.rect(j-1, -12, j, rectval[j],
                col = colval[j])
    circos.text(j-0.5, floor(ylim[2])-2,
                labels = textval[j], cex = 0.70, niceFacing = TRUE)
 }
}
################################################################################
#                                 LINKS                                        #
################################################################################
links$common1 <- links$common
links$common2 <- links$common
for (k in 1:nrow(links))
{
 links$common1[k] <- links$common1[k] / ds[as.numeric(substr(links$s1[k],2,2)),
                                           colnames(ds)==links$GO[k]]
 links$common2[k] <- links$common2[k] / ds[as.numeric(substr(links$s2[k],2,2)),
                                           colnames(ds)==links$GO[k]]
}
# Rescale number of genes in common from 0 to 0.5 (max width):
links$common1 <- links$common1 / 2
links$common2 <- links$common2 / 2
# Colour of links (taken from the previous track):
colv <- col2rgb(colval)

for (i in 1:nrow(links))
{
 posGO <- which(GOs == links$GO[i])
 pos1 <- as.numeric(substr(links$s1[i], 2, 2))
 pos2 <- as.numeric(substr(links$s2[i], 2, 2))
 color <- rgb(colv[1,posGO], colv[2,posGO], colv[3,posGO],
              150, maxColorValue = 255)

 circos.link(pos1,
             c(posGO - 0.5 - links$common1[i], posGO - 0.5 + links$common1[i]),
             pos2,
             c(posGO - 0.5 - links$common2[i], posGO - 0.5 + links$common2[i]),
             col = color, border = "grey", rou = 0.725)
}
################################################################################
#                           FINAL ANNOTATIONS                                  #
################################################################################
circos.rect(-0.05, 0.45, 1, 2.0, xpd = TRUE,
            sector.index = 5, track.index = 1)
circos.text(0.46, 2.3, labels = "Links width",
            sector.index = 5, track.index = 1, cex = 0.75)
circos.text(0.49, 1.6, labels = "0%", cex = 0.70, niceFacing = TRUE,
            sector.index = 5, track.index = 1)
circos.text(0.07, 1.6, labels = "100%", cex = 0.70, niceFacing = TRUE,
            sector.index = 5, track.index = 1)
circos.text(0.87, 1.6, labels = "100%", cex = 0.70, niceFacing = TRUE,
            sector.index = 5, track.index = 1)
circos.text(0.46, 0.8, labels = "of terms in cluster",
            sector.index = 5, track.index = 1, cex = 0.75)
circos.segments(x0 = 0.20, y0 = 1.55, x1 = 0.40 , y1 = 1.55,
                sector.index = 5, track.index = 1)
circos.segments(x0 = 0.56, y0 = 1.55, x1 = 0.76, y1 = 1.55,
                sector.index = 5, track.index = 1)
# Clean up all modified parameters:
circos.clear()
