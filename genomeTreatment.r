# This script was used to prepare genomeLength and genomeLength2000
rawt <- read.table("GCF_000297895.1.assembly.txt", header = TRUE, fill = TRUE)

out <- data.frame(Scaffold = row.names(rawt),
                     Length = rawt$Sequence.Length)

# Remove last row:
out <- out[-nrow(out), ]

write.table(out, file = "genomeLength.txt", row.names = FALSE,
          quote = FALSE, sep = "\t", col.names = FALSE)

################################################################################
ds <- read.table("genomeLength.txt", col.names = c("Scaffold","Length"))

ordre <- order(ds$Scaffold, decreasing = TRUE)
ds <- ds[ordre, ]

totLength <- sum(ds$Length)
ds$Cum <- cumsum(ds$Length)
ds$Prop <- round(((ds$Cum / totLength) * 100), 2)

write.table(ds[1:2000,1:2], file = "genomeLength2000.txt", row.names = FALSE,
          quote = FALSE, sep = "\t", col.names = FALSE)
