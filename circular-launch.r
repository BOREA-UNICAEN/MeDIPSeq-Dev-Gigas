feature <- "CDS"

GOn <- "BP"
source("circular-plot.r")
title(GOn, line = -2)
dev.copy2pdf(device = x11, file = paste(feature,"circle",GOn,".pdf", sep = ""))
GOn <- "MF"
source("circular-plot.r")
title(GOn, line = -2)
dev.copy2pdf(device = x11, file = paste(feature,"circle",GOn,".pdf", sep = ""))
GOn <- "CC"
source("circular-plot.r")
title(GOn, line = -2)
dev.copy2pdf(device = x11, file = paste(feature,"circle",GOn,".pdf", sep = ""))
