"LoadForVariants" <- function(feature, func = "co.var")
{
"co.var" <- function(x) ( sd(x)/mean(x) )

message(paste("Reading source data for ",feature,"...", sep = ""))
# Read the gene universe:
univ <- read.table(paste("counts",feature,"filtered.txt", sep = ""),
                   sep = " ", dec = ".", header = TRUE)
colnames(univ) <- substr(colnames(univ), 1, nchar(colnames(univ))-4)

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

univ <- na.omit(apply(univ, 1, FUN = func))

return(univ)
}
