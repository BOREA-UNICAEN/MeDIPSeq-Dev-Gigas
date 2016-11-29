#gen <- readRDS("matrixReloaded.data")
dmr <- readRDS("allDMR.data")
stages <- c("2-8 cells", "Intermediate", "Spat")
A <- B <- C <- NULL
Agenes <- Bgenes <- Cgenes <- NULL
################################################################################
message("Scanning A...")
subdmr <- subset(dmr, Stage == "2-8 cells")
# For each of the three DMR comparisons (subset):
for (j in 1:nrow(subdmr))
{
 # Identify which scaffold and extract the map table:
 smap <- gen[[ which(names(gen) == subdmr$Item[j]) ]]$map
 # DMR start and stop:
 dmrstart <- subdmr$Start[j]
 dmrstop  <- subdmr$Stop[j]
 # For each feature having a start/stop included in the DMR:
 readstart <- dmrstart < smap$Stop
 readstart <- min(which(readstart))
 readstop <- dmrstop > smap$Start
 readstop <- max(which(readstop))
 if (readstop < readstart)
 { A <- c(A, "OUT") } else {
  # Check if the length of the selection is zero:
  if (!is.infinite(readstart) & !is.infinite(readstop)) {
   # Save that feature on a vector:
   A <- c(A, as.character(smap[readstart:readstop, ]$Feature))
   # Extract the gene IDs for those positions:
   Agenes <- c(Agenes, smap[readstart:readstop, ]$Annot)
  }
 }
}

greppata <- grep("LTR|LINE|DNA|RC", A)
A[greppata] <- "TE"
greppata <- grep("similarity|Repeat", A)
A[greppata] <- "REP"

A <- A[ A != "mRNA" ]
################################################################################
message("Scanning B...")
subdmr <- subset(dmr, Stage == "Intermediate")
# For each of the three DMR comparisons (subset):
for (j in 1:nrow(subdmr))
{
 # Identify which scaffold and extract the map table:
 smap <- gen[[ which(names(gen) == subdmr$Item[j]) ]]$map
 # DMR start and stop:
 dmrstart <- subdmr$Start[j]
 dmrstop  <- subdmr$Stop[j]
 # For each feature having a start/stop included in the DMR:
 readstart <- dmrstart < smap$Stop
 readstart <- min(which(readstart))
 readstop <- dmrstop > smap$Start
 readstop <- max(which(readstop))
 if (readstop < readstart)
 { B <- c(B, "OUT") } else {
  # Check if the length of the selection is zero:
  if (!is.infinite(readstart) & !is.infinite(readstop)) {
   # Save that feature on a vector:
   B <- c(B, as.character(smap[readstart:readstop, ]$Feature))
   # Extract the gene IDs for those positions:
   Bgenes <- c(Bgenes, smap[readstart:readstop, ]$Annot)
  }
 }
}
greppata <- grep("LTR|LINE|DNA|RC", B)
B[greppata] <- "TE"
greppata <- grep("similarity|Repeat", B)
B[greppata] <- "REP"

B <- B[ B != "mRNA" ]
################################################################################
message("Scanning C...")
subdmr <- subset(dmr, Stage == "Spat")
# For each of the three DMR comparisons (subset):
for (j in 1:nrow(subdmr))
{
 # Identify which scaffold and extract the map table:
 smap <- gen[[ which(names(gen) == subdmr$Item[j]) ]]$map
 # DMR start and stop:
 dmrstart <- subdmr$Start[j]
 dmrstop  <- subdmr$Stop[j]
 # For each feature having a start/stop included in the DMR:
 readstart <- dmrstart < smap$Stop
 readstart <- min(which(readstart))
 readstop <- dmrstop > smap$Start
 readstop <- max(which(readstop))
 if (readstop < readstart)
 { C <- c(C, "OUT") } else {
  # Check if the length of the selection is zero:
  if (!is.infinite(readstart) & !is.infinite(readstop)) {
   # Save that feature on a vector:
   C <- c(C, as.character(smap[readstart:readstop, ]$Feature))
   # Extract the gene IDs for those positions:
   Cgenes <- c(Cgenes, smap[readstart:readstop, ]$Annot)
  }
 }
}

greppata <- grep("LTR|LINE|DNA|RC", C)
C[greppata] <- "TE"
greppata <- grep("similarity|Repeat", C)
C[greppata] <- "REP"

C <- C[ C != "mRNA" ]

################################################################################
full <- c(34987154,181064503,25674748,21879155,24463465,233068965)
full <- matrix(nrow = 6, ncol = 1, data = full)
rownames(full) <- c("CDS","intrn","promoter","REP","TE","OUT")
full <- full/sum(full)*100
################################################################################
ds <- data.frame(Counts = c(table(A), table(B), table(C), full),
                 Feature = c(names(table(A)), names(table(B)),
                             names(table(C)), rownames(full)),
                 Comp = c(rep("S",6),rep("I",5),rep("M",6),rep("genome",6)) )

# This re-ordering was changed to take into account that S-I-M are not in the
# alphabetical order:
#ds <- rbind(ds, c(0, "promoter", "I"))
#ds <- ds[order(ds$Comp), ]
ds <- rbind(ds[1:7, ], c(0, "OUT", "I"), ds[8:nrow(ds), ])

ds$Counts <- as.numeric(ds$Counts)

# By comparison, percentage of its maximum:
ds[1:6  , 1] <- with(ds[1:6  , ], Counts/sum(Counts)*100)
ds[7:12 , 1] <- with(ds[7:12 , ], Counts/sum(Counts)*100)
ds[13:18, 1] <- with(ds[13:18, ], Counts/sum(Counts)*100)

# This line added because levels are reordered alphabetically by default:
#levels(ds$Comp) <- unique(as.character(ds$Comp))
ds <- round(xtabs(ds), 3)

write.table(ds, "DMRannotation.csv", sep = ";", dec = ",")
