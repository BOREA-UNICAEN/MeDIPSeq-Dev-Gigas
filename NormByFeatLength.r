############## Normalisation vs feature length #############
rm(list=ls(all=TRUE))
##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/Guillaume_SourceFiles") 
getwd()


#### Load data 
int <- read.table("Cgigas_v9_intron.gff", header = FALSE, sep = "")
cds <- read.table("Cgigas_v9_exon.gff", header = FALSE, sep = "")
gen <- read.table("Cgigas_v9_gene.gff", header = FALSE, sep = "")

# format INT len data data table
INTlen <- data.frame (int[,9], int[,5]-int[,4])
colnames (INTlen) <- c ("GeneID", "length")
INTlen$GeneID <-  substr(as.character(INTlen$GeneID), 8, nchar(as.character(INTlen$GeneID))-1)
# sums INT length for each gene
INTlen <- aggregate(INTlen$length, by = list(INTlen$GeneID), FUN=sum)
colnames (INTlen) <- c ("GeneID", "INTlength")

#does the same for CDS
CDSlen <- data.frame (int[,9], int[,5]-int[,4])
colnames (CDSlen) <- c ("GeneID", "length")
CDSlen$GeneID <-  substr(as.character(CDSlen$GeneID), 8, nchar(as.character(CDSlen$GeneID))-1)
CDSlen <- aggregate(CDSlen$length, by = list(INTlen$GeneID), FUN=sum)
colnames (CDSlen) <- c ("GeneID", "CDSlength")

#does the same for GEN
GENlen <- data.frame (int[,9], int[,5]-int[,4])
colnames (GENlen) <- c ("GeneID", "length")
GENlen$GeneID <-  substr(as.character(GENlen$GeneID), 8, nchar(as.character(GENlen$GeneID))-1)
GENlen <- aggregate(GENlen$length, by = list(INTlen$GeneID), FUN=sum)
colnames (GENlen) <- c ("GeneID", "GENlength")