################################# Bar Plot Reads In Features Distribution Rev Fig1a ###########################
rm(list=ls(all=TRUE))
##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/sommes de comptes normalises") 
getwd()
########Chargement des dpnnées 
message("loading data...")
library(ggplot2)
#Distribution of Counts in Features
counts <- read.table("TableOfSumsFinal.txt", header=T, sep="\t")
counts <- counts[1:5,]
colnames(counts)=c("Feature", colnames(counts[,2:22]))
counts[,2:22] <- cumsum(counts[,2:22])
counts <- data.frame (counts[,1], 
						(counts[,2]+counts[,4])/2,
						(counts[,5]+counts[,6]+counts[,7])/3,
						(counts[,8]+counts[,9])/2,
						(counts[,10]+counts[,11])/2,
						(counts[,12]+counts[,13]+counts[,14])/3,
						(counts[,15]+counts[,16])/2,
						(counts[,17]+counts[,18]+counts[,19])/3,
						(counts[,20]+counts[,21]+counts[,22])/3
                      )
colnames(counts)=c("Feature", "ovo", "C", "Mor", "Bla", "Gas", "Tro", "Dla", "Spa")					  
Ovo <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$ovo[1:4]),ymax=c(counts$ovo[1:5]))
C <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$C[1:4]),ymax=c(counts$C[1:5]))
Mor <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Mor[1:4]),ymax=c(counts$Mor[1:5]))
Bla <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Bla[1:4]),ymax=c(counts$Bla[1:5]))
Gas <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Gas[1:4]),ymax=c(counts$Gas[1:5]))
Tro <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Tro[1:4]),ymax=c(counts$Tro[1:5]))
Dla <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Dla[1:4]),ymax=c(counts$Dla[1:5]))
Spa <- data.frame(Feature=c("CDS", "INT", "PRO", "REP", "TE"), ymin=c(0, counts$Spa[1:4]),ymax=c(counts$Spa[1:5]))

################################################################################
stages <- c("Ovocyte","2-8 cells","Morula","Blastula","Gastrula",
            "Trocophore","D-larva","Spat")
# PREPARE THE COLOUR VECTORS FOR THE PLOT:
oricols <- c("aquamarine3","cadetblue3","chartreuse3","darkorchid3","deeppink3")
rgbs <- col2rgb(oricols)
hsvs <- rgb2hsv(rgbs)
valueseq <- seq(0.4, 0.9, length.out = 8)

HSV <- data.frame(h = rep(hsvs[1,], each = 8),
                  s = rep(hsvs[2,], each = 8),
                  v = rep(valueseq, 5))
HSV <- t(HSV)

colvec <- character(40)
for(i in 1:40) colvec[i] <- hsv(HSV[1,i], HSV[2,i], HSV[3,i])
colmat <- matrix(nrow = 8, ncol = 5, byrow = FALSE, data = colvec)
rownames(colmat) <- stages
colnames(colmat) <- c("CDS","INT","PRO","REP","TE")

p<- ggplot(Ovo, aes(fill=feature,ylab = "Proportion of reads mapped in features", ymax=ymax, ymin=ymin, xmax=0.90, xmin=0.1)) +
     geom_rect(colour="black", fill=colmat[1,]) +
	 geom_rect(data=C, xmax=1.90, xmin=1.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[2,]) +
	 geom_rect(data=Mor, xmax=2.9, xmin=2.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[3,]) +
	 geom_rect(data=Bla, xmax=3.9, xmin=3.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[4,]) +
	 geom_rect(data=Gas, xmax=4.9, xmin=4.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[5,]) +
	 geom_rect(data=Tro, xmax=5.9, xmin=5.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[6,]) +
	 geom_rect(data=Dla, xmax=6.9, xmin=6.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[7,]) +
	 geom_rect(data=Spa, xmax=7.9, xmin=7.1, aes(ymax=ymax, ymin=ymin),colour="black", fill=colmat[8,]) +
     #coord_polar(theta="y") +
     xlim(c(0, 8)) + 
	 #theme(panel.grid=element_blank()) +
	 #theme(panel.background=element_blank()) +
     #theme(axis.text=element_text("feature")) +
     #theme(axis.ticks=element_blank()) +
	 theme(aspect.ratio=4) +
     labs(title="Distribution of mapped reads in features", xlab="Stage", ylab="proportion of reads mapped in features")
save.image(p, "revfig1a", png)


############# Peasron's chi square test for distribution
#chargement du fichier contenant le nb de counts
counts <- read.table("TableOfSumsFinalRPKM.txt", header=T, sep="\t")

counts <- data.frame (counts[,1], 
						(counts[,2]+counts[,4])/2,
						(counts[,5]+counts[,6]+counts[,7])/3,
						(counts[,8]+counts[,9])/2,
						(counts[,10]+counts[,11])/2,
						(counts[,12]+counts[,13]+counts[,14])/3,
						(counts[,15]+counts[,16])/2,
						(counts[,17]+counts[,18]+counts[,19])/3,
						(counts[,20]+counts[,21]+counts[,22])/3
                      )
colnames(counts)=c("Feature", "ovo", "C", "Mor", "Bla", "Gas", "Tro", "Dla", "Spa")					  

Xsq <- chisq.test(
					counts[,2:9], correct = TRUE,
					p = rep(1/length(counts[,8:9]), length(counts[,8:9])), rescale.p = TRUE
				)

### Kolmogorov Smirnof for pairwise comparisons between de stages
ks <- ks.test (counts$ovo, counts$C)

#### WAYYYYYY TO COMPLICATED AND UGLY !!!! WHAT SHOUL DBE USED IS ######
## counts <- read.table("TableOfSumsFinal.txt", header=T, sep="\t")
## colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts))-1)
## bystage <- aggregate(t(counts), by = list(rownames(t(counts))), FUN = mean)
## rownames(bystage) <- bystage[ ,1]
##counts <- t(bystage[ ,-1])
### Reorder:
## ordre <- c("ovo","X28cell","morula","blastula","gastrula","troc","dlarv","spat")
## counts <- counts[ ,ordre]