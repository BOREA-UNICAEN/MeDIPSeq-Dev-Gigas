suppressPackageStartupMessages( library(topGO) )

stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")
GOs <- c("BP","MF","CC")

message(paste("Reading source data ..."))
# Read the gene universe:
univ <- read.table("DMRs/expressions.csv", header = T, sep = ";", dec = ",")[,1]
# The file GOmappings.txt has been previously parsed by a PERL script:
readmap <- readMappings("GOmappings.txt")
# Read the vector that assigns some genes to being methylated:
kmc <- readRDS("K-means.data")$cluster

# Prepare a vector that links each gene to its cluster number, or to zero:
universe <- vector("numeric", length = length(univ))
names(universe) <- univ
universe[names(universe) %in% names(kmc)] <- 1

# Prepare a selection function:
"ourSelectionFunction" <- function(allScore) {
   return(allScore == 1)
}

################################################################################
message("Finding most enriched GOs...")

for (i in 1:3)                         # For each gene ontology...
{
  # Build topGOdata object:
  GOdata <- new("topGOdata", description = "Methylated vs universe",
                ontology = GOs[i], allGenes = universe,
                geneSel = ourSelectionFunction,
                annot = annFUN.gene2GO, gene2GO = mappings)

  # Run tests:
  #testKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

  # Table of most enriched GOs for that stage and that ontology:
  enrich <- GenTable(GOdata, Fvalue = testF, topNodes = 10)
  print(enrich)
}
