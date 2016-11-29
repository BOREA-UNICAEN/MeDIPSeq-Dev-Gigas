source("C:/Users/Samuele/Documents/Dropbox/Antroposeine/LIM modelling/phenagen.r")
source("LoadForVariants.r")
################################################################################
features <- c("PRO","CDS","INT")
vary <- read.table("variants.txt", header = TRUE)
#selection <- readRDS("significantGenes.data")

# Filtering away single-varying gene IDs:
vary <- vary[vary$variants != 1, ]
vary <- vary[vary$variants < 40, ]

par(mfrow = c(1, 3), mar = c(4, 4, 3, 0.5))

# For each feature:
for (i in 1:length(features))
{
 ds <- LoadForVariants(features[i])
 # Keep only the GeneIDs that are present in the filtered variants table:
 ds <- ds[names(ds) %in% vary$GeneID]

 # Keep only the GeneIDs present in the selection vector:
 #ds <- ds[names(ds) %in% selection]

 model <- gamen(x = vary$variants[vary$GeneID %in% names(ds)],
                y = ds,
                xlab = "Number of gene variants",
                ylab = "Variability of methylation",
                plot = TRUE, fam = gaussian, ylim = c(0, 2))
 title(main = features[i])
}
################################################################################
nonmethylated <- !(vary$GeneID %in% rownames(meth))
plusvary <- vary[!nonmethylated, ]
moinvary <- vary[nonmethylated, ]


dat <- data.frame(Condition = c(rep("Meth",nrow(plusvary)),
                                rep("Non-Meth",nrow(moinvary))),
                  Variants = c(plusvary$variants, moinvary$variants))
ggplot(dat, aes(x=Condition, y=Variants, fill=Condition)) +
 geom_boxplot(outlier.shape = TRUE) + ylim(0, 5)
