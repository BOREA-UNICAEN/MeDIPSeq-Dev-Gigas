suppressPackageStartupMessages( library(topGO) ) # car fait chier

stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")
GOs <- c("BP","MF","CC")
feature <- "CDS"

# Read the gene universe:
message(paste("Reading source data for ",feature,"...", sep = ""))
univ <- read.table(paste("counts",feature,"filtered.txt", sep = ""),
                   sep = " ", dec = ".", header = TRUE)
# The file GOmappings.txt has been previously parsed by a PERL script:
mappings <- readMappings("GOmappings.txt")
# Selection vector of interesting GOs:
selection <- as.vector(read.csv("selectedGOs.csv")[,1])

colnames(univ) <- substr(colnames(univ), 1, nchar(colnames(univ))-4)   #########

# ANOVA by row (factor is "development stage"):
message("Running ANOVA against development stages...")
testRes <- numeric(nrow(univ))
for (i in 1:nrow(univ))
{
 values <- as.numeric(univ[i, ])
 test <- aov(values ~ as.factor(colnames(univ)))
 testRes[i] <- summary(test)[[1]]$"Pr(>F)"[1]
}
# Exclude all genes rows for which ANOVA gave non-significant results:
univ <- univ[testRes < 0.05, ]

genesToUse <- rownames(univ)
# Reads the other universe:
univ <- read.csv("LCRcountsCDS.csv", header = TRUE)
# Keep in the LCR table only those genes corresponding to the ANOVA selection:
univ <- univ[univ$X %in% genesToUse, ]
# Do the always useful mean by stages:
rownames(univ) <- univ$X; univ <- univ[ ,-1]
colnames(univ) <- substr(colnames(univ), 1, nchar(colnames(univ))-1)   #########
message("Calculating means by development stages...")
bystage <- aggregate(t(univ), by = list(rownames(t(univ))), FUN = mean)
rownames(bystage) <- bystage[ ,1]
univ <- t(bystage[ ,-1])
univ <- univ[ ,stages]

# Prepare a selection function:
"ourSelectionFunction" <- function(allScore) return (allScore)

################################################################################
# Empty variables to store results:
methres <- numeric(5)
selmat <- matrix(0, nrow = 14, ncol = 8)
scolist <- NULL

message("Finding most enriched GOs...")
i = 1

for (j in 1:length(stages))           # For each development stage:
 {
  # Select one development stage:
  vect <- univ[ ,stages[j]]

  # Build topGOdata object:
  invisible (capture.output(
   GOdata <- new("topGOdata", description = stages[j],
                 ontology = GOs[i],
                 allGenes = vect, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)    ))
  # Run tests:
  invisible (capture.output(
   testKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
   #testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  ))
  # Table of FIVE most enriched GOs for that stage and that ontology:
  #enrich <- GenTable(GOdata, KSvalue = testKS, topNodes = 5)
  # List of all genes with those FIVE most enriched GOs:
  #enrichedgenes <- genesInTerm(GOdata, enrich$GO.ID)

  #vec <- score(testKS, whichGO = selection)
  #names(vec) <- Term(names(vec))
  #selmat[ ,j] <- vec

  scol <- scoresInTerm(GOdata, selection)
  if (is.list(scolist)) {
   for (k in 1:length(scolist)) scolist[[k]] <- c(scolist[[k]], scol[[k]])
  } else { scolist <- scol }

}
saveRDS(scolist, "selectedGOmeth.data")
################################################################################
# Empty variables to store results:
selmat <- matrix(0, nrow = 14, ncol = 8)
colnames(selmat) <- stages

i = 1

for (j in 1:length(stages))           # For each development stage:
 {
  # Select one development stage:
  vect <- univ[ ,stages[j]]

  # Build topGOdata object:
  invisible (capture.output(
   GOdata <- new("topGOdata", description = stages[j],
                 ontology = GOs[i],
                 allGenes = vect, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)    ))

  # List of genes with the GOs identified:
  whichgenes <- genesInTerm(GOdata, names(scolist))

  for (k in 1:length(whichgenes)) selmat[k,j] <- mean ( vect[whichgenes[[k]]] )
}
rownames(selmat) <- names(whichgenes)

saveRDS(selmat, "selectedGOlcr.data")
