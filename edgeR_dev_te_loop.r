#This script processes in edgeR the counts files output from htseq-count.
library("edgeR")

#Get the list of counts files to process.
setwd("counts_dev/te")

listfiles <- list.files()
#Read the files and store in a DGEList object.
data = readDGE(listfiles)

########## Filter out weakly expressed and non-informative features ###########

print("Filtering data...")
#Create a logical vector noint where the rows of data containing the __no_feature and other uninformative data are flagged as TRUE in the vector.
noint = rownames(data) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
#Calculate counts per million.
cpmd = cpm(data)
#Keep if counts per million is greater than one for at least two of the samples and the value of noint is not FALSE (i.e. = TRUE).
#This will remove the __no_feature data at the end of the counts files.
#In at least two of the samples because we have a minimum of 2 datasets per condition.
keep = rowSums(cpmd > 1) >=2 & !noint
#Now reset counts to contain only the values we want to keep.
data = data[keep,]
#Reset the library sizes:
data$samples$lib.size = colSums(data$counts)

########## Analyze MDS plots ###########

#Associate groups to the sample names using the file targets.txt.
targets <- readTargets(file="../targets.txt")
# Remove ovo3 from targets:
targets <- targets[targets$Sample != "ovo3" , ]
data <- DGEList(counts = data, group=targets$Stage)
#Normalize data to reduce competition effects.
print("Normalizing...")
data = calcNormFactors(data)

########## Do pairwise comparisons ###########

print("Estimating dispersions...")
#Estimate dispersions.
data = estimateCommonDisp(data)
data = estimateTagwiseDisp(data)

#Detect significantly different counts between groups using the exact binomial test.
# Sam: this part was modified to perform comparisons in order.
if (yespairs == TRUE) {
print("Doing pair-wise comparisons...")

stages <- c("ovocyte","28cell","morula","blastula","gastrula","trocophore",
            "dlarvae","spat")
allcombs <- combn(stages, 2)

extest = list()
k = 1

for (i in 1:ncol(allcombs))
{
 extest[[k]] <- exactTest(data, pair = c(allcombs[1,i], allcombs[2,i]))
 k <- k + 1
}

}