suppressPackageStartupMessages( library(topGO) )

stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")
GOs <- c("BP","MF","CC")
feature <- "GEN"

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

# Completely inelegant way of calculating means by development stage:
message("Calculating means by development stages...")
bystage <- aggregate(t(univ), by = list(rownames(t(univ))), FUN = mean)
rownames(bystage) <- bystage[ ,1]
univ <- t(bystage[ ,-1])
univ <- univ[ ,stages]

# Prepare a selection function:
"ourSelectionFunction" <- function(allScore) return (allScore)

#"ourSelectionFunction" <- function(allScore) {
#   return(allScore > 0)
#}

################################################################################
# Empty table to store results:
result <- data.frame(stage = rep(stages, each = 5),
                     BP = "", BP_Score = 0, MF = "", MF_Score = 0,
                     CC = "", CC_Score = 0,
                     stringsAsFactors = FALSE)
methres <- numeric(5)
selmat <- matrix(0, nrow = 14, ncol = 8)
scolist <- NULL

message("Finding most enriched GOs...")
for (i in 1:3)                         # For each gene ontology...
{
 for (j in 1:length(stages))           # ... and for each development stage.
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
  enrich <- GenTable(GOdata, KSvalue = testKS, topNodes = 5)
  # List of all genes with those FIVE most enriched GOs:
  enrichedgenes <- genesInTerm(GOdata, enrich$GO.ID)
  # Mean score of those genes:
  for (k in 1:5) {
      methvalues <- geneScore(GOdata, whichGenes = unlist(enrichedgenes[[k]]))
      methres[k] <- mean(methvalues)
  }

  #vec <- score(testKS, whichGO = selection)
  #names(vec) <- Term(names(vec))
  #selmat[ ,j] <- vec

  scol <- scoresInTerm(GOdata, selection)
  if (is.list(scolist)) {
   for (k in 1:length(scolist)) scolist[[k]] <- c(scolist[[k]], scol[[k]])
  } else { scolist <- scol }

  # Extract enrichment results into a data.frame:
  result[result$stage == stages[j] , colnames(result) == GOs[i]] <- enrich$GO.ID
  result[result$stage == stages[j] , which(colnames(result) == GOs[i])+1] <- methres
 }
}

saveRDS(result, paste("enrich-methyl-",feature,".data", sep = ""))
ds <- result

################################################################################
#       Translates GOs to their complete names for a nice presentation:        #
################################################################################
totrows <- c(rep("BP", length(unique(result$BP))),
             rep("MF", length(unique(result$MF))),
             rep("CC", length(unique(result$CC))))
totcodes <- c(unique(result$BP), unique(result$MF), unique(result$CC))
totdefs <- Term(totcodes)

nicetable <- data.frame(Feature = feature,
                        Ontology = totrows,
                        Code     = totcodes,
                        Definition = totdefs,
                        row.names = NULL)

################################################################################
# Create another matrix with the number of genes in common for the selected GO #
################################################################################
for (j in 1:3)
{
 message(paste("Analysing genes in common for",GOs[j],"ontology...", sep = " "))
 liston <- split(ds[ ,colnames(ds) == GOs[j]], ds$stage)
 combs <- combn(names(liston), 2)

 # For each combination of stages, retrieve a list (common) of the intersecting
 # GOs, then build two vectors of stage names (with length determined by the
 # number of intersecting GOs).
 common <- list(NULL)
 firstname <- secondname <- character()
 links <- list(NULL)

 for (i in 1:ncol(combs))
 {
  common[[i]] <- intersect(liston[[combs[1,i]]], liston[[combs[2,i]]])
  firstname <- c(firstname, rep(combs[1,i], times = length(common[[i]])))
  secondname <- c(secondname, rep(combs[2,i], times = length(common[[i]])))
 }

 # The following data structure will be used by the circular plot script:
 links <- data.frame(s1 = firstname, s2 = secondname,
                     GO = unlist(common), common = 0,
                     stringsAsFactors = FALSE)

 for (i in 1:nrow(links))
 {
  message(paste("Iteration:",i,"of",nrow(links)))
  vect1 <- univ[ ,colnames(univ) == links$s1[i]]
  vect2 <- univ[ ,colnames(univ) == links$s2[i]]

  invisible (capture.output(
  GOdata1 <- new("topGOdata", description = "First stage",
                 ontology = GOs[j],
                 allGenes = vect1, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)  ))
  invisible (capture.output(
  GOdata2 <- new("topGOdata", description = "Second stage",
                 ontology = GOs[j],
                 allGenes = vect2, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)  ))
  # Extract the list of genes included in the GO term:
  genes1 <- unlist(genesInTerm(GOdata1, links$GO[i]))
  genes2 <- unlist(genesInTerm(GOdata2, links$GO[i]))
  # Number of common genes is equal to length of intersection:
  links$common[i] <- length(intersect(genes1, genes2))
 }
 saveRDS(links, paste("genes_common_",GOs[j],".",feature,".data", sep = ""))
}