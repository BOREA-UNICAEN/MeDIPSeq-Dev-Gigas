extract <- data.frame(Scaffold = "",
                      Annot = "",
                      Distance = 0,
                      Expr = 0,
                      GeneLength = 0)

# Extract only mRNA features:
for (i in 1:length(gen))
{
 submrna <- subset(gen[[names(gen)[i]]]$map, Feature == "mRNA")
 # Delete rows for which expression is undefined:
 submrna <- na.omit(submrna)

 if (!is.null(submrna$DistA28c)) distanze<-submrna$DistA28c else distanze<-rep(NA, nrow(submrna))
 submrna <- cbind(submrna, distanze)

 # Extract expression levels:
 scaffExpr <- expr$logA[expr$GeneID %in% submrna$Annot]
 names(scaffExpr) <- expr$GeneID[expr$GeneID %in% submrna$Annot]
 #scaffExpr <- expr[expr$GeneID %in% submrna$Annot, ]

 # Filter out genes for which we have no expression levels (very few):
 submrna <- submrna[submrna$Annot %in% names(scaffExpr), ]

 extract <- rbind(extract, data.frame(Scaffold = rep(names(gen)[i],nrow(submrna)),
                                      Annot = submrna$Annot,
                                      Distance = submrna$distanze,
                                      Expr = scaffExpr,
                                      GeneLength = abs(submrna$Start - submrna$Stop)))
}

# Clean out first row:
extract <- extract[-1, ]
extract <- do.call(data.frame,lapply(extract, function(x) replace(x, is.infinite(x),NA)))
extract <- na.omit(extract)

extract <- subset(extract, Expr != 0)

# Plots:
#xs <- extract$Expr[extract$Distance < 100000]
#ys <- extract$Expr[extract$Distance > 100000]
#facteur <- c(rep("lower", length(xs)), rep("higher", length(ys)))
#fisher.test(x = c(xs,ys), y = facteur)
extract$Distance <- extract$Distance / 1000

smoothScatter(extract$Distance, extract$Expr,
              nbin = 200,
              colramp = colorRampPalette(c("white","blue2","red")),
              axes = FALSE, frame.plot = TRUE, #ylim = c(-9.6, +9.6),
              xlab = "", ylab = "Gene expression (logFC)")
axis(side = 1, las = 1)
axis(side = 2, las = 1)
mtext("Distance to DMR (kbp)", side = 1, line = 1.75, cex = 0.63)

#plot(extract$Distance, extract$Expr, col = "grey80", pch = 19,
#     main = "",
#     xlab = "", ylab = "Gene expression logFC",
#     ylim = c(-10, +10), axes = FALSE, frame.plot = TRUE)
#axis(side = 1, las = 1)
#axis(side = 2, las = 1)
#y1 <- par("usr")[3]; y2 <- par("usr")[4]
#par(new = TRUE)
#bplot <- boxplot(extract$Expr ~ quantiles, outline = FALSE, axes = FALSE,
#        ylim = c(y1, y2), at = apply(limits, 1, mean))
