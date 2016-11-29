suppressPackageStartupMessages( library(topGO) )

GOs <- c("BP","MF","CC")

message(paste("Reading source data ..."))
# Read the gene universe:
univ <- read.table("DMRs/expressions.csv", header = T, sep = ";", dec = ",")[,1]
# The file GOmappings.txt has been previously parsed by a PERL script:
mappings <- readMappings("GOmappings.txt")

# Read the vector that assigns some genes to clusters (not all of them):
kmc <- readRDS("K-means.data")$cluster

# Prepare a vector that links each gene to its cluster number, or to zero:
universe <- vector("numeric", length = length(univ))
names(universe) <- univ
universe[names(universe) %in% names(kmc)] <- kmc

# Prepare a selection function:
"ourSelectionFunction" <- function(allScore) {
   return(allScore == j)
}

################################################################################
# Empty table to store results:
result <- data.frame(cluster = 1:5,
                     BP = as.numeric(0), MF = as.numeric(0), CC = as.numeric(0),
                     stringsAsFactors = FALSE)
resgo <- list(BP = "", MF = "", CC = "")

message("Finding most enriched GOs...")
for (i in 1:3)                         # For each gene ontology...
{
 resclus <- list(C1 = "", C2 = "", C3 = "", C4 = "", C5 = "")
 for (j in 1:5)                 # ... and for each cluster.
 {
  # Build topGOdata object:
  invisible(capture.output(
   GOdata <- new("topGOdata", description = paste("Cluster",j),
                 ontology = GOs[i], allGenes = universe,
                 geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings) ))
  # Run tests:
  invisible (capture.output(
   #testKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
   testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  ))

  # Table of most enriched GOs for that stage and that ontology:
  enrich <- GenTable(GOdata, Fvalue = testF, topNodes = 200)
  # Select only those with Fvalue < 0.05:
  enrich <- enrich[as.numeric(enrich$Fvalue) < 0.05, ]

  result[result$cluster == j,colnames(result) == GOs[i]] <- nrow(enrich)

  resclus[[j]] <- enrich$GO.ID
  # Extract results into a data.frame:
  #result[result$cluster == j, colnames(result) == GOs[i]] <- enrich$GO.ID
 }
 resgo[[i]] <- resclus
}

#result[ ,2:4] <- as.numeric(result[ ,2:4])
saveRDS(result, "enrich-methyl-CLUSTERS.data")
################################################################################
#       Translates GOs to their complete names for a nice presentation:        #
################################################################################
#totrows <- c(rep("BP", length(unique(result$BP))),
#             rep("MF", length(unique(result$MF))),
#             rep("CC", length(unique(result$CC))))
#totcodes <- c(unique(result$BP), unique(result$MF), unique(result$CC))
#totdefs <- Term(totcodes)

#nicetable <- data.frame(Feature = feature,
#                        Ontology = totrows,
#                        Code     = totcodes,
#                        Definition = totdefs,
#                        row.names = NULL)
################################################################################
# Create another matrix with the number of genes in common for the selected GO #
################################################################################
totlinks <- data.frame(s1 = "", s2 = "",
                     GO = "", common = 0,
                     stringsAsFactors = FALSE)

for (j in 1:3)
{
 message(paste("Analysing terms in common for",GOs[j],"ontology...", sep = " "))

 liston <- resgo[[j]]
 combs <- combn(names(liston), 2)

# For each combination of clusters, retrieve a 'liston' of the intersecting GOs,
# then build two vectors of stage names (with length determined by the number of
# intersecting GOs).
 common <- list(NULL)
 firstname <- secondname <- character()
 links <- list(NULL)

 for (i in 1:ncol(combs))
 {
  common[[i]] <- intersect(liston[[combs[1,i]]], liston[[combs[2,i]]])
  firstname <- c(firstname, rep(combs[1,i], times = length(common[[i]])))
  secondname <- c(secondname, rep(combs[2,i], times = length(common[[i]])))
 }

 links <- data.frame(s1 = combs[1, ], s2 = combs[2, ],
                     GO = GOs[j], common = 0,
                     stringsAsFactors = FALSE)

 # Number of common genes is equal to length of intersection:
 for (i in 1:nrow(links)) links$common[i] <- length(common[[i]])

 totlinks <- rbind(totlinks, links)
}
totlinks <- totlinks[-1, ]
saveRDS(totlinks, "genes_common_clusters.data")
