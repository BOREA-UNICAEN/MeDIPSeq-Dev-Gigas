ds <- read.table("DMRannotation.csv", sep = ";", dec = ",", header = TRUE,
                 row.names = 1)
ds <- as.matrix(ds)

colvec <- c("aquamarine3","cadetblue3","chartreuse3","darkorchid3","deeppink3",
            "white")
#par(mar = c(2, 3.5, 2.5, 4.5), xpd = TRUE, mgp = c(2.5, 1, 0))
par(mgp = c(2.5, 0.75, 0))
barx <- barplot(ds, space = c(0.2, 0.2, 0.2, 0.5), #beside = FALSE,
          axes = FALSE, main = "",
          col = colvec,
          ylab = "DMR annotation (%)",
          ylim = c(0, 100))
axis(side = 2, las = 1)
mtext(c("N = 1043","N = 14","N = 2230",""), at = barx, side = 1,
      cex = 0.6, line = 1.8)
################################################################################

