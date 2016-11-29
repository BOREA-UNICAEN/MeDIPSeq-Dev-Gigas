ds <- readRDS("wholeGenome.data")
shortlist <- read.table("genomeLength2000.txt", header = FALSE)
colnames(shortlist) <- c("scaffold","length")
################################################################################
#   Select only those scaffolds that were selected into the 2000 shortlist     #
################################################################################
message("Selecting scaffolds in the shortlist...")
# Empty logical vector to store results:
IsIn <- logical(nrow(ds))

for (i in 1:nrow(ds))
 {
  # Is the first column of the i-th row present in the selection file? If so,
  # isin[i] will be TRUE.
  IsIn[i] <- ds$scaffold[i] %in% shortlist$scaffold
 }

# Keep only rows corresponding to TRUE:
ds <- ds[IsIn, ]

message(paste(length(unique(ds$scaffold)))," scaffolds kept")
################################################################################
message("Processing scaffolds...")
allscaffs <- unique(ds$scaffold)

# Empty list:
liston <- vector("list", length(allscaffs))
names(liston) <- allscaffs

for (i in 1:length(allscaffs))
{
 # First prepare the map table to be inserted:
 select <- ds[ds$scaffold == allscaffs[i], ]
 select <- select[order(select$start), 2:ncol(select)]

 # Clean up the annotation string:
 select$annot <- gsub(";", "", select$annot)
 select$annot <- gsub("Parent=", "", select$annot)
 select$annot <- gsub("ID=", "", select$annot)

 # Change column names to be as awesome as possible:
 colnames(select) <- c("Feature","Start","Stop","Strand","Annot")

 # Add expression value next to each identified gene ID:
 #select$Expr <- expdata$value[match(select$Annot, expdata$GeneID)]

 # Select which is the transcription start:
 select$TranStart <- ifelse(select$Strand == "+", select$Start, select$Stop)

 # Finally, insert it at its correct position in the list:
 liston[[i]]$map <- select
}
################################################################################
#      Preparing the second table with differential methylation regions        #
################################################################################
allscaffs <- as.vector(allscaffs)
dmrfile <- c("drA.txt","drB.txt","drC.txt")
stages <- c("A28c","BInt","CSpa")
allDMR <- NULL

for (f in 1:length(dmrfile))
{
 message(paste("Processing DMRs for stage: ", stages[f], sep = ""))
 dr <- read.table(paste("Dev experiment/diffRepsOutput/", dmrfile[f], sep = ""),
                  header = TRUE)
 # Filter only those those DMR with a strong p-value:
 dr <- subset(dr, pval < 1e-7)
 # List of scaffolds to be processed for this particular diffReps comparison:
 dmrscaffs <- unique(dr$Chrom)

 for (s in 1:length(dmrscaffs))
 {
  # Select one scaffold at a time:
  select <- subset(dr, Chrom == dmrscaffs[s])

  # Prepare a nice table for the output:
  newframe <- data.frame(Item = select$Chrom,   # Leave this in case of debug
                         Start = select$Start,
                         Stop = select$End,
                         Length = select$Length,
                         logFC = select$log2FC,
                         Stage = stages[f])
  # Find the position in the list that hosts the same scaffold:
  listpos <- which(allscaffs %in% dmrscaffs[s])

  # Continue only if the scaffold is found in the genome:
  if (length(listpos) != 0)
  {
   # Remove DMRs for which we have infinite logFC values:
   newframe[!is.finite(newframe$logFC), ] <- NA
   newframe <- na.omit(newframe)

   allDMR <- rbind(allDMR, newframe)
   # If its position in the list is not already occupied by a dmr table, then
   # create it; else, append it at the end:
   if (is.null(liston[[listpos]]$dmr)) { liston[[listpos]]$dmr <- newframe
      } else {
      liston[[listpos]]$dmr <- rbind(liston[[listpos]]$dmr, newframe) }
  }
 }
}
################################################################################
#      Calculate distance to nearest DMR for each genetic feature              #
################################################################################
message("Calculating distances to nearest DMR centres...")

for (i in 1:length(liston))
{
 distances <- numeric(nrow(liston[[i]]$map))
 relDMRfc  <- numeric(nrow(liston[[i]]$map))
 # If this scaffold has a dmr table, continue. Otherwise, move to next scaffold.
 if (!is.null(liston[[i]]$dmr))
 {
  for (s in 1:length(stages))
  {
   # Extract from the DMR table only the records relative to the iterated stage:
   DMR <- liston[[i]]$dmr[liston[[i]]$dmr$Stage == stages[s], ]
   # The position of a DMR is assumed as its mean point:
   DMRpos <- (DMR$Start + DMR$Stop) / 2
   DMRpos <- round(DMRpos)
   # Then, cycle through the features in the scaffold:
   for (j in 1:nrow(liston[[i]]$map))
   {
    # Vector of all possible distances of that feature from all DMRs:
    allPossDistances <- abs(liston[[i]]$map$TranStart[j] - DMRpos)
    # Retrieve the minimum of those distances:
    distances[j] <- min(allPossDistances)
    # Save the logFC of the DMR that was just selected as the nearest one:
    # (the MIN function was added to avoid solving my programming problems)
    relDMRfc[j] <- min( DMR$logFC[allPossDistances == min(allPossDistances)] )
   }   # END FOR EACH FEATURE
   # Append the new data to the map table:
   liston[[i]]$map <- cbind(liston[[i]]$map, distances, relDMRfc)
   # Put the correct names to the new columns:
   colnames(liston[[i]]$map)[ncol(liston[[i]]$map)-1] <- paste("Dist",stages[s],sep="")
   colnames(liston[[i]]$map)[ncol(liston[[i]]$map)]   <- paste("Fold",stages[s],sep="")
  }    # END FOR EACH STAGE
 }     # END IF
}      # END FOR EACH SCAFFOLD
