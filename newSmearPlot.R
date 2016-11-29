thisdir <- getwd()

yespairs <- TRUE
source("counts_dev/edgeR_dev_cds_loop.r")

# Vector of stages order:
stages <- c("ovocyte","28cell","morula","blastula","gastrula","trocophore",
            "dlarvae","spat")
stagenames <- c("Ovocyte","2-8 cells","Morula","Blastula","Gastrula",
                "Trocophore","D-larva","Spat")

split.screen(c(9, 9), erase = TRUE)
seqrow <- seq(1, 81, 9)

for(i in 1:length(extest))
{
   sig <- decideTestsDGE(extest[[i]], p = 0.01)
   sigtags <- rownames(data)[as.logical(sig)]

   pos1 <- which(stages %in% extest[[i]]$comparison[1])
   pos2 <- which(stages %in% extest[[i]]$comparison[2])

   position <- seqrow[pos1+1] + pos2
   screen(position)
   par(mar = c(0.5, 0.5, 0.5, 0.5))

   plotSmear(extest[[i]], de.tags = sigtags,
             allCol = "black", deCol = "#64C8A6",
             axes = FALSE, ylim = c(-8.0, +8.0))
   axis(side = 2, cex.axis = 0.6, las = 1, at = c(-8, 0, +8), labels = FALSE)
   axis(side = 1, at = c(2, 6, 12), labels = FALSE, tcl = -0.3)
   if (position %in% (seqrow+1)) {
       mtext(c("-8","0","+8"), side = 2, at = c(-8, 0, +8), cex = 0.55, las = 1,
             line = 0.6) }
}

################################################################################
setwd(thisdir)
source("counts_dev/edgeR_dev_te_loop.r")

for(i in 1:length(extest))
{
   sig <- decideTestsDGE(extest[[i]], p = 0.01)
   sigtags <- rownames(data)[as.logical(sig)]

   pos1 <- which(stages %in% extest[[i]]$comparison[1])
   pos2 <- which(stages %in% extest[[i]]$comparison[2])

   position <- seqrow[pos2+1] + pos1
   screen(position)
   par(mar = c(0.5, 0.5, 0.5, 0.5))

   plotSmear(extest[[i]], de.tags = sigtags,
             allCol = "black", deCol = "#E61284",
             axes = FALSE, xlim = c(2, 12), ylim = c(-8.0, +8.0))
   axis(side = 2, cex.axis = 0.6, las = 1, at = c(-8, 0, +8), labels = FALSE)
   axis(side = 1, at = c(2, 6, 12), labels = FALSE, tcl = -0.3)
   #if (position %in% (seqrow+1)) {
   #   mtext(c("-8","0","+8"), side = 2, at = c(-8, 0, +8), cex = 0.55, las = 1,
   #         line = 0.6) }
}
################################################################################

# Stage names:
for (i in 1:8)
{
 screen(i+1)
 par(mar = c(0.5, 0.5, 0.5, 0.0))
 plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, frame.plot = FALSE)
 text(0.5, 0.5, labels = stagenames[i], cex = 0.85)
}

for (i in 1:8)
{
 screen(seqrow[i+1])
 par(mar = c(0.5, 0.5, 0.5, 0.5))
 plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, frame.plot = FALSE)
 text(0.47, 0.5, labels = stagenames[i], cex = 0.85)
}

# Legend panel: (MIGHT NEED TO EDIT IT IN POST-PRODUCTION)
screen(1)
par(mar = c(2,2,0.2,0.2))
#par(mar = c(0.2, 0.2, 0.2, 0.2))
plot(NULL, xlim = c(1, 12), ylim = c(-8,+8), axes = FALSE, frame.plot = TRUE)
axis(side = 1, cex.axis = 0.6, at = c(2, 6, 12), labels = FALSE)
mtext(c("2","6","12"), side = 1, at = c(2, 6, 12), cex = 0.65, las = 1,
      line = 0.35)
axis(side = 2, cex.axis = 0.6, las = 1, at = c(-8, 0, +8), labels = FALSE)
mtext(c("-8","0","+8"), side = 2, at = c(-8, 0, +8), cex = 0.65, las = 1,
      line = 0.6)
mtext("log CPM", side = 1, cex = 0.7, line = 1.0)
mtext("log FC", side = 2, cex = 0.7, line = 1.2)

#text(6.5, 0, labels = "Column\nvs\nrow", cex = 0.6)

setwd(thisdir)
close.screen(all = TRUE)