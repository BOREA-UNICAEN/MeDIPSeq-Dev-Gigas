library(xlsx)
gen <- readRDS("matrixReloaded.data")
expr <- read.xlsx2("expressions.xlsx", 1, header = TRUE, colClasses =
                  c("character", rep("numeric",15)))
# Select expression columns:
expr <- subset(expr, select = c(1, 14, 15, 16))

################################################################################
par(mfrow = c(3, 1), mar = c(2, 4, 0.5, 2))
ylimits <- c(-15, +15)

source("DMRtreatmentA.r")
mtext("A", side = 4, line = 0.3, at = 10, las = 1, cex = 1.2)
source("DMRtreatmentB.r")
mtext("B", side = 4, line = 0.3, at = 15, las = 1, cex = 1.2)
par(mar = c(4, 4, 0.5, 2))
source("DMRtreatmentC.r")
mtext("C", side = 4, line = 0.3, at = 15, las = 1, cex = 1.2)
