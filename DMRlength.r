"replot" <- function(col = "grey80", ylab = "")
{
  hist(subd$Length, freq = FALSE,
      xlim = c(1000, 6000), ylim = c(0, 10e-04),
      main = "", xlab = "", ylab = "",
      col = col, axes = FALSE)
  mtext("DMR length (bp)", side = 1, line = 1.75, cex = 0.63)
  mtext(ylab,              side = 2, line = 2.8,  cex = 0.63)
  box()
  axis(side = 1, las = 1)
  axis(side = 2, las = 1)
}

ds <- readRDS("allDMR.data")

#par(mfrow = c(3, 1), mar = c(3.1, 4.2, 0.2, 0.2), mgp = c(3.3, 0.8, 0))

A <- subd <- subset(ds, Stage == "2-8 cells")
replot("grey30", ylab = "Frequency")
text(5995, 9.8e-04, labels = "S", cex = 1.4)

B <- subd <- subset(ds, Stage == "Intermediate")
replot("grey50")
text(5995, 9.8e-04, labels = "I", cex = 1.4)

C <- subd <- subset(ds, Stage == "Spat")
replot("grey80")
text(6000, 9.8e-04, labels = "M", cex = 1.4)

#print(ks.test(A$Length, B$Length))
#print(ks.test(B$Length, C$Length))
#print(ks.test(A$Length, C$Length))
