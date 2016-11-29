message(paste("*** Feature is ",feature," ***", sep = ""))
ds <- read.csv(paste("LCRcounts",feature,".csv", sep = ""),
               header = TRUE, row.names = 1)

# Substitute stage names for versions good for plotting:
stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")
"switchnumber" <- function (x) switch(x, "ovo"   = 1,
                                      "X28cell"  = 2,
                                      "morula"   = 3,
                                      "blastula" = 4,
                                      "gastrula" = 5,
                                      "troc"     = 6,
                                      "dlarv"    = 7,
                                      "spat"     = 8, x)
# THE USUAL TREATMENT FOR CDS METHYLATION DATA #################################
colnames(ds) <- substr(colnames(ds), 1, nchar(colnames(ds))-1)                 #
message("Running ANOVA against development stages...")                         #
testRes <- numeric(nrow(ds))                                                   #
for (i in 1:nrow(ds))                                                          #
{                                                                              #
 values <- as.numeric(ds[i, ])                                                 #
 test <- aov(values ~ as.factor(colnames(ds)))                                 #
 testRes[i] <- summary(test)[[1]]$"Pr(>F)"[1]                                  #
}                                                                              #
ds <- ds[testRes < 0.0001, ]                                                   #
message(paste(dim(ds)[1], " genes kept", sep = ""))                            #
message("Calculating means by development stages...")                          #
bystage <- aggregate(t(ds), by = list(rownames(t(ds))), FUN = mean)            #
rownames(bystage) <- bystage[ ,1]                                              #
ds <- t(bystage[ ,-1])                                                         #
# Reorder matrix for stages:                                                   #
ds <- ds[ ,stages]                                                             #
# Update stage names:                                                          #
for (i in 1:ncol(ds)) colnames(ds)[i] <- switchnumber(colnames(ds)[i])         #
# ##############################################################################
# Hierarchical clustering:
hc <- hclust(dist(ds), method = "average")
# Re-order lines according to clustering:
ds <- ds[hc$order, ]

tri <- ds
rownames(tri) <- 1:dim(ds)[1]
################################################################################
# Matrix to plot an awesome bottom gradient:
gradmat <- matrix(nrow = nrow(tri), ncol = ncol(tri)*8, byrow = FALSE,
                  data = rep(1:64, each = nrow(tri)))

library(plot3D)
colpal <- colorRampPalette(c("blue2","blue2","grey75",
                             "firebrick1","firebrick1"))(256)
grapal <- colorRampPalette(c("grey25","grey95"))(256)
#zfacet <- tri[-1, -1] + tri[-1, -ncol(tri)] + tri[-nrow(tri), -1] + tri[-nrow(tri), -ncol(tri)]
#facetcol <- cut(zfacet, 256)
persp3D(z = tri, phi = phi, theta = theta, resfac = 10, bty = "b2",
      xlab = xlabel, ylab = "Development", zlab = "Methylation",
      main = feature,
      col = colpal, border = NA, colkey = FALSE, cex.lab = 1.5)

image3D(y = seq(0, 0.995, length.out = ncol(gradmat)), z = min(tri),
        col = grapal, colvar = gradmat,
        add = TRUE, colkey = FALSE)