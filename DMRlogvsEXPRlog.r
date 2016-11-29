#library(xlsx)
#gen <- readRDS("matrixReloaded.data")
#expr <- read.xlsx2("expressions.xlsx", 1, header = TRUE, colClasses =
#                  c("character", rep("numeric",15)))
# Select expression columns:
#expr <- subset(expr, select = c(1, 14, 15, 16))

# Graphic parameters:
#par(mfcol = c(3, 1), mar = c(3.5, 3.8, 0.5, 1.3), mgp = c(2.5, 0.9, 0))
colour <- c("grey65","firebrick1","blue2")

extract <- data.frame(#Scaffold = "",
                      Annot = "",
                      DMRfc = 0,
                      GENfc = 0)
# Empty list to store contingency tables:
counts <- list()

# Extract only mRNA features:
for (i in 1:length(gen))
{
 submrna <- subset(gen[[names(gen)[i]]]$map, Feature == "mRNA")
 # Delete rows for which expression is undefined:
 submrna <- na.omit(submrna)

 if (!is.null(submrna$FoldA28c)) logfc<-submrna$FoldA28c else logfc<-rep(NA, nrow(submrna))
 submrna <- cbind(submrna, logfc)

 # Filter out genes not in the selection list:
 #submrna <- submrna[submrna$Annot %in% selection, ]

 # Extract expression levels:
 scaffExpr <- expr$logA[expr$GeneID %in% submrna$Annot]
 names(scaffExpr) <- expr$GeneID[expr$GeneID %in% submrna$Annot]
 #scaffExpr <- expr[expr$GeneID %in% submrna$Annot, ]

 # Filter out genes for which we have no expression levels (very few):
 submrna <- submrna[submrna$Annot %in% names(scaffExpr), ]


 if (nrow(submrna) != 0) {
 extract <- rbind(extract, data.frame(#Scaffold = rep(names(gen)[i],nrow(submrna)),
                                      Annot = submrna$Annot,
                                      DMRfc = submrna$logfc,
                                      GENfc = scaffExpr)) }
}

extract <- extract[-1, ]
extract <- do.call(data.frame,lapply(extract, function(x) replace(x, is.infinite(x),NA)))
work <- na.omit(extract)

colvec <- character(nrow(work))
colvec[work$GENfc == 0] <- colour[1]
colvec[work$GENfc > 0]  <- colour[2]
colvec[work$GENfc < 0]  <- colour[3]
rgbvec <- col2rgb(colvec)
colvec <- rgb(rgbvec[1,],rgbvec[2,],rgbvec[3,], 95, maxColorValue = 255)

plot(work$DMRfc, work$GENfc, pch = 19, cex = 0.8, axes = FALSE, frame.plot = T,
     col = colvec,
               xlim = c(-9.3, +5), ylim = c(-10, +16),
               xlab = "", ylab = "Gene expression (logFC)")
axis(side = 1, las = 1)
axis(side = 2, las = 1)
abline(v = 0, lty = 2, col = "grey65")
text(4.9, 15.3, labels = "S", cex = 1.4)
#mtext("A", side = 4, line = 0.2, at = 16, las = 1, cex = 1.2)
mtext("DMR methylation (logFC)", side = 1, line = 1.75, cex = 0.63)

# Add the text for the previous plot:
text(4.9, 52.5, labels = "S", cex = 1.4, xpd = NA)

# Once plot is done, calculate contingency table:
work$DMRfc[work$DMRfc > 0] <- +1
work$DMRfc[work$DMRfc < 0] <- -1
work$GENfc[work$GENfc > 0] <- +1
work$GENfc[work$GENfc < 0] <- -1
counts$A <- with(work, table(DMRfc, GENfc))
################################################################################
extract <- data.frame(#Scaffold = "",
                      Annot = "",
                      DMRfc = 0,
                      GENfc = 0)

# Extract only mRNA features:
for (i in 1:length(gen))
{
 submrna <- subset(gen[[names(gen)[i]]]$map, Feature == "mRNA")
 # Delete rows for which expression is undefined:
 submrna <- na.omit(submrna)

 if (!is.null(submrna$FoldBInt)) logfc<-submrna$FoldBInt else logfc<-rep(NA, nrow(submrna))
 submrna <- cbind(submrna, logfc)

 # Filter out genes not in the selection list:
 #submrna <- submrna[submrna$Annot %in% selection, ]

 # Extract expression levels:
 scaffExpr <- expr$logB[expr$GeneID %in% submrna$Annot]
 names(scaffExpr) <- expr$GeneID[expr$GeneID %in% submrna$Annot]
 #scaffExpr <- expr[expr$GeneID %in% submrna$Annot, ]

 # Filter out genes for which we have no expression levels (very few):
 submrna <- submrna[submrna$Annot %in% names(scaffExpr), ]



 if (nrow(submrna) != 0) {
 extract <- rbind(extract, data.frame(#Scaffold = rep(names(gen)[i],nrow(submrna)),
                                      Annot = submrna$Annot,
                                      DMRfc = submrna$logfc,
                                      GENfc = scaffExpr)) }
}

extract <- extract[-1, ]
extract <- do.call(data.frame,lapply(extract, function(x) replace(x, is.infinite(x),NA)))
work <- na.omit(extract)

colvec <- character(nrow(work))
colvec[work$GENfc == 0] <- colour[1]
colvec[work$GENfc > 0]  <- colour[2]
colvec[work$GENfc < 0]  <- colour[3]
rgbvec <- col2rgb(colvec)
colvec <- rgb(rgbvec[1,],rgbvec[2,],rgbvec[3,], 95, maxColorValue = 255)

plot(work$DMRfc, work$GENfc, pch = 19, cex = 0.8, axes = FALSE, frame.plot = T,
     col = colvec,
               xlim = c(-9.3, +5), ylim = c(-10, +16),
               xlab = "", ylab = "")
axis(side = 1, las = 1)
axis(side = 2, las = 1)
abline(v = 0, lty = 2, col = "grey65")
text(4.9, 15.3, labels = "I", cex = 1.4)
#mtext("B", side = 4, line = 0.2, at = 16, las = 1, cex = 1.2)
mtext("DMR methylation (logFC)", side = 1, line = 1.75, cex = 0.63)

# Add the text for the previous plot:
text(4.9, 52.5, labels = "I", cex = 1.4, xpd = NA)

# Once plot is done, calculate contingency table:
work$DMRfc[work$DMRfc > 0] <- +1
work$DMRfc[work$DMRfc < 0] <- -1
work$GENfc[work$GENfc > 0] <- +1
work$GENfc[work$GENfc < 0] <- -1
counts$B <- with(work, table(DMRfc, GENfc))
################################################################################
extract <- data.frame(#Scaffold = "",
                      Annot = "",
                      DMRfc = 0,
                      GENfc = 0)

# Extract only mRNA features:
for (i in 1:length(gen))
{
 submrna <- subset(gen[[names(gen)[i]]]$map, Feature == "mRNA")
 # Delete rows for which expression is undefined:
 submrna <- na.omit(submrna)

 if (!is.null(submrna$FoldCSpa)) logfc<-submrna$FoldCSpa else logfc<-rep(NA, nrow(submrna))
 submrna <- cbind(submrna, logfc)

 # Filter out genes not in the selection list:
 #submrna <- submrna[submrna$Annot %in% selection, ]

 # Extract expression levels:
 scaffExpr <- expr$logC[expr$GeneID %in% submrna$Annot]
 names(scaffExpr) <- expr$GeneID[expr$GeneID %in% submrna$Annot]
 #scaffExpr <- expr[expr$GeneID %in% submrna$Annot, ]

 # Filter out genes for which we have no expression levels (very few):
 submrna <- submrna[submrna$Annot %in% names(scaffExpr), ]

 if (nrow(submrna) != 0) {
 extract <- rbind(extract, data.frame(#Scaffold = rep(names(gen)[i],nrow(submrna)),
                                      Annot = submrna$Annot,
                                      DMRfc = submrna$logfc,
                                      GENfc = scaffExpr)) }
}

extract <- extract[-1, ]
extract <- do.call(data.frame,lapply(extract, function(x) replace(x, is.infinite(x),NA)))
work <- na.omit(extract)

colvec <- character(nrow(work))
colvec[work$GENfc == 0] <- colour[1]
colvec[work$GENfc > 0]  <- colour[2]
colvec[work$GENfc < 0]  <- colour[3]
rgbvec <- col2rgb(colvec)
colvec <- rgb(rgbvec[1,],rgbvec[2,],rgbvec[3,], 95, maxColorValue = 255)

#par(mar = c(3.5, 3.8, 0.5, 1.3))
plot(work$DMRfc, work$GENfc, pch = 19, cex = 0.8, axes = FALSE, frame.plot = T,
     col = colvec,
               xlim = c(-9.3, +5), ylim = c(-10, +16),
               xlab = "", ylab = "")
axis(side = 1, las = 1)
axis(side = 2, las = 1)
abline(v = 0, lty = 2, col = "grey65")
text(4.9, 15.3, labels = "M", cex = 1.4)
#mtext("C", side = 4, line = 0.2, at = 16, las = 1, cex = 1.5)
mtext("DMR methylation (logFC)", side = 1, line = 1.75, cex = 0.63)

# Add the text for the previous plot:
text(4.9, 52.5, labels = "M", cex = 1.4, xpd = NA)

# Once plot is done, calculate contingency table:
work$DMRfc[work$DMRfc > 0] <- +1
work$DMRfc[work$DMRfc < 0] <- -1
work$GENfc[work$GENfc > 0] <- +1
work$GENfc[work$GENfc < 0] <- -1
counts$C <- with(work, table(DMRfc, GENfc))

#chisq.test(counts$A[ ,-2])
#chisq.test(counts$B[ ,-2])
#chisq.test(counts$C[ ,-2])