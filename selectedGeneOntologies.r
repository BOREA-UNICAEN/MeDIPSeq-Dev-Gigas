suppressPackageStartupMessages( library(GO.db) )
source("Plotranges.R")
std <- function(x) sd(x)/sqrt(length(x))
matr <- readRDS("selectedGOtest.data")
ds <-   readRDS("selectedGOmeth.data")
lcr <-  readRDS("selectedGOlcr.data")

matr <- matr[-1, ]
ordermanual <- c(8,1,9,11,14,3,10,12,7,13,2,4,6,5)
matr <- matr[rev(ordermanual), ]

rownames(matr)[rownames(matr) ==
  "regulation of transcription, DNA-templated"] <- "regulation of transcription"
################################################################################
# We need to prepare a data frame with min, max, and central values:
rang <- data.frame(Name = character(length(ds)),
                   Min = numeric(length(ds)),
                   Central = numeric(length(ds)),
                   Max = numeric(length(ds)),
                   N = numeric(length(ds)),
                   stringsAsFactors = FALSE)
for (i in 1:length(ds))
{
 rang$Name[i] <- names(ds)[i]
 rang$Central[i] <- mean(ds[[i]])
 rang$Min[i] <- mean(ds[[i]]) - std(ds[[i]])
 rang$Max[i] <- mean(ds[[i]]) + std(ds[[i]])
 rang$N[i]   <- length(ds[[i]])
}
rang$Name <- Term(rang$Name)
rang <- rang[rev(ordermanual), ]
#rang$Min[rang$Min < 0] <- 0
################################################################################
#ll <- layout(matrix(ncol = 4, c(1,1,2,3), byrow = TRUE))
# creates a color palette
my_palette <- colorRampPalette(c("lightskyblue4", "lightskyblue3", "grey90"))(3)
col_breaks = c(seq(0   ,0.05, length = 1),
               seq(0.05,0.1 , length = 1),
               seq(0.1 ,1   , length = 2))

par(mar = c(4.8, 0.8, 0.5, 13.6))
image(z = t(matr),
      col = my_palette,
      breaks = col_breaks,
      axes = FALSE)
# Names of stages:
text(seq(0, 1, length.out = ncol(matr)), par("usr")[3]-0.022,
     labels = colnames(matr),
     srt = 70, adj = 1, xpd = TRUE, cex = 1)
# Names of ontologies:
text(1.08, seq(0, 1, length.out = nrow(matr)),
     labels = rownames(matr), adj = 0, xpd = TRUE, cex = 1.0)

# Add the p-values on each cell:
for (i in 1:nrow(matr)) {
 for (j in 1:ncol(matr)) {
  xposit <- seq(0, 1, length.out = ncol(matr))[j]
  yposit <- seq(0, 1, length.out = nrow(matr))[i]
  text(xposit, yposit, labels = round(matr[i,j], 2))
 }
}

rowseq <- seq(0, 1, length.out = ncol(matr))
rowGridlines <- rowseq[2 : (length(rowseq))] - (rowseq[2] / 2)
abline(v = rowGridlines, col = "grey", lwd = 1.2)
colseq <- seq(0, 1, length.out = nrow(matr))
colGridlines <- colseq[2 : (length(colseq))] - (colseq[2] / 2)
abline(h = colGridlines, col = "grey", lwd = 1.2)
################################################################################
rownames(lcr) <- Term(rownames(lcr))
lcr <- lcr[rev(ordermanual), ]
# Create a color palette:
my_palette <- colorRampPalette(c("blue2", "grey85", "firebrick3"))(299)
col_breaks = c(seq(min(lcr),-0.031, length = 100),
               seq(0.051   , 0.1  , length = 100),
               seq(0.031   , 1    , length = 100))

par(mar = c(4.8, 0.5, 0.5, 0.0))
image(z = t(lcr),
      col = my_palette,
      #breaks = col_breaks,
      axes = FALSE)
text(seq(0, 1, length.out = ncol(lcr)), par("usr")[3]-0.022,
     labels = colnames(matr),
     srt = 70, adj = 1, xpd = TRUE, cex = 1)
#text(1.08, seq(0, 1, length.out = nrow(lcr)),
#     labels = rownames(lcr), adj = 0, xpd = TRUE, cex = 0.75)

rowseq <- seq(0, 1, length.out = ncol(lcr))
rowGridlines <- rowseq[2 : (length(rowseq))] - (rowseq[2] / 2)
abline(v = rowGridlines, col = "grey", lwd = 1.2)
colseq <- seq(0, 1, length.out = nrow(lcr))
colGridlines <- colseq[2 : (length(colseq))] - (colseq[2] / 2)
abline(h = colGridlines, col = "grey", lwd = 1.2)
################################################################################
par(mar = c(4.8, 2.3, 0.5, 0.7))
Plotranges(min = rang$Min, max = rang$Max, value = rang$Central,
           labels = rang$N, xlab = "Methylation (counts)")
