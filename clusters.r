stagenames <- c("Ovocyte","2-8 cells","Morula","Blastula","Gastrula",
                "Trocophore","D-larva","Spat")
kmc <- readRDS("K-means.data")
ds <- read.csv("LCRcountsCDS.csv", header = TRUE)
rownames(ds) <- ds[ ,1]; ds <- ds[ ,-1]
colnames(ds) <- substr(colnames(ds), 1, nchar(colnames(ds))-1)

# ANOVA by row (factor is "development stage"):
message("Running ANOVA against development stages...")
testRes <- numeric(nrow(ds))
for (i in 1:nrow(ds))
{
 values <- as.numeric(ds[i, ])
 test <- aov(values ~ as.factor(colnames(ds)))
 testRes[i] <- summary(test)[[1]]$"Pr(>F)"[1]
}

# Exclude all genes rows for which ANOVA gave non-significant results:
ds <- ds[testRes < 0.0001, ]
message(paste(dim(ds)[1], " genes kept", sep = ""))

# Completely inelegant way of calculating means by development stage:
bystage <- aggregate(t(ds), by = list(rownames(t(ds))), FUN = mean)
rownames(bystage) <- bystage[ ,1]
bystage <- t(bystage[ ,-1])
# Reorder the table by columns:
bystage <- bystage[ ,c(5, 8, 4, 1, 3, 7, 2, 6)]

# Some filtering because reasons:
greppata <- bystage[ ,"morula"] > -2;   bystage <- bystage[greppata, ]
greppata <- bystage[ ,"blastula"] > -2; bystage <- bystage[greppata, ]
greppata <- bystage[ ,"dlarv"] > -2;    bystage <- bystage[greppata, ]

# Colour palette:
col_fun <- colorRamp(c("blue2","grey70","firebrick1"))
# Quick function to rescale each value vector to 0-1 range:
"reScale" <- function(x, from, to) (x - min(x)) / max(x - min(x)) * (to - from) + from
# Generate the colour vector:
colmat <- reScale(bystage, 0, 1)

# Plot:
#par(mfrow = c(1, 5))
for (i in 1:5)
{
 GenesInClus <- names( kmc$cluster[kmc$cluster == i] )
 subm <- bystage[rownames(bystage) %in% GenesInClus, ]
 subc <- colmat [rownames(colmat)  %in% GenesInClus, ]
 plot(NULL, xlim = c(1,8), ylim = c(-3.0, +3.2),
      axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
 axis(side = 1, las = 2, at = 1:8, labels = NA)
 axis(side = 2, las = 1)
 text(4.5, 3.1, labels = paste("Cluster",i))
 for (j in 1:nrow(subm))
 {
  colm <- col_fun(subc[j,1:8])
  colm <- rgb(colm[,1], colm[,2], colm[,3], 255, maxColorValue = 255)
  lines(1:8, subm[j,1:8], col = "grey80")
  points(1:8, subm[j,1:8], col = colm, pch = 19)
 }
}
################################################################################
#GenesThatChange <- list(names(kmc$cluster)[kmc$cluster == 1],
#                        names(kmc$cluster)[kmc$cluster == 2],
#                        names(kmc$cluster)[kmc$cluster == 3],
#                        names(kmc$cluster)[kmc$cluster == 4],
#                        names(kmc$cluster)[kmc$cluster == 5])

#down28   <- names(kmc$cluster)[kmc$cluster == 1 | kmc$cluster == 4]
#downSpat <- names(kmc$cluster)[kmc$cluster == 5]
#upSpat   <- names(kmc$cluster)[kmc$cluster == 1 | kmc$cluster == 3]
#GenesThatChange <- list("down28" = down28, "downSpat" = downSpat,
#                        "upSpat" = upSpat)
#saveRDS(GenesThatChange, file = "GenesThatChange.data")