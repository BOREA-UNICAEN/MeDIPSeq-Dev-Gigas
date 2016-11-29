# This script iteratively reads all *.bed files in the working directory. For
# each one, it creates a logical vector indicating which rows are present in
# a selection file (TRUE), or absent (FALSE). It then extracts only the rows
# with TRUE, and saves the resulting filtered table in the output file.

filelist <- list.files()
# Keep only files with "bed" in their name:
filelist <- filelist[grep("bed", filelist)]
# Selection file:
genome <- read.table("genomeLength2000.txt", sep = "\t", header = FALSE)[ ,1]

for (f in 1:length(filelist))
{
 # For each file in the list, read it:
 message(paste("Reading ",filelist[f],"...", sep = ""))
 bed <- read.table(filelist[f], sep = "\t", header = FALSE)
 message(paste(nrow(bed), "rows"))
 
 # Empty logical vector to store results:
 isin <- logical(nrow(bed))

 message("Going through the data...")
 for (i in 1:nrow(bed))
 {
  # Is the first column of the i-th row present in the selection file? If so,
  # isin[i] will be TRUE.
  isin[i] <- bed[i,1] %in% genome
 }

 message("Extracting and saving file...")
 # Keep only rows corresponding to TRUE:
 bed <- bed[isin, ]
 # Saves output, adding an X before the original name:
 write.table(bed, sep = "\t",
             col.names = FALSE, row.names = FALSE, quote = FALSE,
             file = paste("X",filelist[f], sep = ""))
 message("*****")
}

