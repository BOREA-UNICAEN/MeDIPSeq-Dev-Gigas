library(topGO)
source("DMRannotation.r")
mappings <- readMappings("C:/Users/Samuele/Documents/Dropbox/Huitres/Dev experiment/GOmappings.txt")

# Clean up the vectors:
greppata <- grep("CGI", Agenes)
Agenes <- Agenes[greppata]
greppata <- grep("CGI", Bgenes)
Bgenes <- Bgenes[greppata]
greppata <- grep("CGI", Agenes)
Cgenes <- Cgenes[greppata]

Agenes <- unique(Agenes);Bgenes <- unique(Bgenes);Cgenes <- unique(Cgenes);

univ <- read.table("expressions.csv", header = T, sep = ";", dec = ",")[,1]
# Prepare a vector that links each gene to its cluster number, or to zero:
universe <- vector("numeric", length = length(univ))
names(universe) <- univ

# Prepare a selection function:
"ourSelectionFunction" <- function(allScore) {
   return(allScore == 1)
}

################################################################################
universe[names(universe) %in% Agenes] <- 1
GOdata <- new("topGOdata", description = "DMR A",
                 ontology = "BP",
                 allGenes = universe, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)
testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
enrich <- GenTable(GOdata, Fvalue = testF, topNodes = 10)
################################################################################
universe[] <- 0
universe[names(universe) %in% Bgenes] <- 1
GOdata <- new("topGOdata", description = "DMR B",
                 ontology = "BP",
                 allGenes = universe, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)
testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
enrich <- GenTable(GOdata, Fvalue = testF, topNodes = 10)
################################################################################
universe[] <- 0
universe[names(universe) %in% Cgenes] <- 1
GOdata <- new("topGOdata", description = "DMR C",
                 ontology = "BP",
                 allGenes = universe, geneSel = ourSelectionFunction,
                 annot = annFUN.gene2GO, gene2GO = mappings)
testF  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
enrich <- GenTable(GOdata, Fvalue = testF, topNodes = 10)







Ag <- summary(as.factor(Agenes))[1:20]
mappos <- which( names(readmap) %in% names(Ag) )
unname(Term(unlist(readmap[mappos])))

Bg <- summary(as.factor(Bgenes))[1:20]
mappos <- which( names(readmap) %in% names(Bg) )
unname(Term(unlist(readmap[mappos])))

Cg <- summary(as.factor(Cgenes))[1:20]
mappos <- which( names(readmap) %in% names(Cg) )
unname(Term(unlist(readmap[mappos])))
