# This scripts requires the 64-bit version of R.
suppressPackageStartupMessages( library(GO.db) )
library(gplots)
stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")

ds <- readRDS("enrich-methyl-GEN.data")

# Substitute stage names for versions good for plotting:
"switchnames" <- function (x) switch(x, "ovo" = "Ovocyte",
                                    "X28cell" = "2-8 cells",
                                     "morula" = "Morula",
                                   "blastula" = "Blastula",
                                   "gastrula" = "Gastrula",
                                       "troc" = "Trocophore",
                                      "dlarv" = "D-larva",
                                       "spat" = "Spat", x)

matr <- data.frame(stage = rep(ds$stage, 3))
matr$GO <- c(ds$BP, ds$MF, ds$CC)
matr$score <- c(ds$BP_Score, ds$MF_Score, ds$CC_Score)

matr <- data.frame(stage = rep(ds$stage, 1))
matr$GO <- ds$BP
matr$score <- ds$BP_Score


# Create contingency table using GO and stage as the two marginal factors:
matr <- xtabs(score ~ GO + stage, data = matr)
# Reorder matrix for stages:
matr <- matr[ ,stages]

for (i in 1:length(colnames(matr))) colnames(matr)[i] <- switchnames(colnames(matr)[i])

# Translate GO into their complete names:
rownames(matr) <- Term(rownames(matr))

################################################################################
# creates a color palette from red to green
my_palette <- colorRampPalette(c("dodgerblue3", "gold1", "firebrick3"),
                               bias = 0.5)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1 ,0  , length = 100),              # for red
               seq(0  ,0.8, length = 100),              # for yellow
               seq(0.8,1  , length = 100))              # for green

heatmap.2(matr,
  trace = "none",         # turns off trace lines inside the heat map
#  tracecol = "black",
  colsep = 0:length(stages)+1, rowsep = 0:nrow(matr),
  srtCol = 60,
  sepcolor = "grey",
  sepwidth = c(0.001, 0.001),
  margins = c(6,15),     # widens margins around plot
  col = my_palette,       # use on color palette defined earlier
#  breaks = col_breaks,    # enable color transition at specified limits
  dendrogram = "row",     # only draw a row dendrogram
  Colv = NA,              # turn off column clustering
  adjCol = 1,
  cexRow = 0.87,
  offsetRow = 0.1, offsetCol = 0.1,
  # Key plot:
  keysize = 1,
  density.info = "none",
  key.xlab = "Methylation score",
  key.title = "",
  key.xtickfun = function() {
        breaks <- parent.frame()$breaks
        return(list(at = parent.frame()$scale01(c(breaks[1],
                                                  breaks[length(breaks)/2],
                                                  breaks[length(breaks)])),
                    labels = c(as.character(breaks[1]),
                               round(breaks[length(breaks)/2], 1),
                               round(breaks[length(breaks)], 1))
              ))
        }
)  # End of heatmap.2 function
