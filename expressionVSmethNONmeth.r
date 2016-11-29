stages <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")

message("Reading expression and methylation data...")
expr <- read.csv("expressionsStagesOnly.csv", header = TRUE)
meth <- read.table("countsCDSfiltered.txt", sep = " ", dec = ".", header = TRUE)

colnames(meth) <- substr(colnames(meth), 1, nchar(colnames(meth))-4)   #########
bystage <- aggregate(t(meth), by = list(rownames(t(meth))), FUN = mean)
rownames(bystage) <- bystage[ ,1]
meth <- t(bystage[ ,-1])

rownames(expr) <- expr[ ,1]
expr <- expr[ ,-1]

nonmethylated <- !(rownames(expr) %in% rownames(meth))
nonmeth <- expr[nonmethylated, ]
meth <- expr[!nonmethylated, ]

meth <- as.vector(as.matrix(meth))
nonmeth <- as.vector(as.matrix(nonmeth))


dat <- data.frame(Condition = c(rep("Non-Meth",length(nonmeth)),
                                rep("Meth",length(meth))),
                  Expression = c(nonmeth, meth))
dat$Condition <- factor(dat$Condition, c("Non-Meth","Meth"))
################################################################################
pal <- colorRampPalette(c("blue2","firebrick1"))

boxplot(dat$Expression ~ dat$Condition,
        outline = FALSE, notch = FALSE,
        las = 1,
        names = c("Non-methylated","Methylated"),
        ylab = "",
        col = pal(2),
        ylim = c(0, 60),
        par(mar = par("mar")),
        pars = list(whisklty = 1, whiskcol = "grey40",
                    staplelty = 1, staplecol = "grey40"))
mtext("Expression", side = 2, line = 2.2, cex = 0.85)
