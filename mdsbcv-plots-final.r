# Grabs this working directory for later use:
thisdir <- getwd()

# Avoid losing time for pair-wide comparisons in the sourced scripts:
yespairs <- FALSE

par(mfrow = c(1, 2), mar = c(4, 4, 3, 0.5), mgp = c(2.5, 0.8, 0))

# Vector of stages order:
stages <- c("ovo","28cell","morula","blastula","gastrula","troc",
            "dlarv","spat")
stagenames <- c("Ovocyte","2-8 cells","Morula","Blastula","Gastrula",
                "Trocophore","D-larva","Spat")
# Vector of colours:
colTE <- c("#37041F","#37041F","#37041F","#710941","#710941","#C81073",
          "#C81073","#C81073","#8E0B52","#8E0B52","#8E0B52","#540730","#540730",
          "#1A020F","#1A020F","#E61284","#E61284","#E61284","#AB0D63","#AB0D63")
colCDS <- c("#1B372D","#1B372D","#1B372D","#38715E","#38715E","#64C8A6",
          "#64C8A6","#64C8A6","#478E76","#478E76","#478E76","#2A5445","#2A5445",
          "#0D1A15","#0D1A15","#72E6BE","#72E6BE","#72E6BE","#55AB8E","#55AB8E")
labvec <- c("2-8","2-8","2-8","BLA","BLA","DLV","DLV","DLV","GAS","GAS","GAS",
            "MOR","MOR","OVO","OVO","SPA","SPA","SPA","TRO","TRO")

# This script returns the data object:
source("counts_dev/edgeR_dev_te_loop.r")

plotMDS(data, col = colTE, labels = labvec, font = 2, method = "bcv",
        xaxt = "n", yaxt = "n", main = "TE")
par(font = 1)
axis(side = 1, las = 1)
axis(side = 2, las = 1)

# Reset working directoy and proceed with the next one:
setwd(thisdir)
source("counts_dev/edgeR_dev_cds_loop.r")

plotMDS(data, col = colCDS, labels = labvec, font = 2, method = "bcv",
        xaxt = "n", yaxt = "n", main = "CDS")
par(font = 1)
axis(side = 1, las = 1)
axis(side = 2, las = 1)

# Reset it again because I need to do trials:
setwd(thisdir)
