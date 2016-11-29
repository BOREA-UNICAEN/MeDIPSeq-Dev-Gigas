############## INT / mRNA expression #############
rm(list=ls(all=TRUE))
##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/Guillaume_SourceFiles") 
getwd()

########Chargement des dpnnées 
message("loading data...")
#variants
variants<- read.table("expdata/variants.txt", header=T, sep="\t")
#mRNA expression par stade (from oyster v9)
mRNAByStage<-read.table("expdata/expdevmerge.txt", header=T, sep="\t")
mRNACV<-read.table("expdata/expdevCV.txt", header=T, sep="\t")
#longueur des gènes pour normalisation
GENlength<-read.table("GENlength.txt", header=T, sep="\t")
#methylation
methCDS<-read.table("countsCDSfiltered.txt", header=T, sep="")
methINT<-read.table("countsINTfiltered.txt", header=T, sep="")
methGB<-read.table("countsGENfiltered.txt", header=T, sep="\t")
FractInCDS<-read.table("FractionInCDSDevAll.txt", header=T, sep="\t")
methCDSvsINT <- read.table("countsnCDS_INT_filtered_devall.txt",header = TRUE, sep="\t")
methPRO<- read.table("countsPROfiltered.txt",header = TRUE, sep="")
#Gene subsets
DevGiga<-read.table("gigaton_dev.txt",header = TRUE, sep="\t")
hox<-read.table("Homeobox_genes.txt",header = TRUE, sep="\t")


######### ANOVA ratioINTCDS vs DEV STAGE ###########
# ANOVA by row (factor is "development stage"):
ds <- methINT
rownames(ds) <- ds[ ,1]; ds <- ds[ ,-1]
colnames(ds) <- substr(colnames(ds), 1, nchar(colnames(ds))-4)
testRes <- numeric(nrow(ds))

message("Running ANOVA against development stages...")

for (i in 1:nrow(ds))
{
 values <- as.numeric(ds[i, ])
 test <- aov(values ~ as.factor(colnames(ds)))
 testRes[i] <- summary(test)[[1]]$"Pr(>F)"[1]
}
methINTDevAll<-data.frame(methINT, testRes)


# Exclude all genes rows for which ANOVA gave non-significant results:
df<-methINTDevAll[!is.na(methINTDevAll$testRes)& methINTDevAll$testRes< 0.01,] 
message(paste(dim(df)[1], " genes kept", sep = ""))



#moyenne par stade colonnes ordonnées
colnames(methINT)
# [1] "GeneID"       "X28cell1pro"  "X28cell2pro"  "X28cell3pro"  "blastula2pro"
 #[6] "blastula3pro" "dlarv1pro"    "dlarv2pro"    "dlarv3pro"    "gastrula1pro"
#[11] "gastrula2pro" "gastrula3pro" "morula1pro"   "morula2pro"   "ovo1pro"     
#[16] "ovo2pro"      "spat1pro"     "spat2pro"     "spat3pro"     "troc1pro"    
#[21] "troc2pro" 

methINTByStage<-data.frame(df[,1],
									(df[,15]+df[,16])/2,#ovo
									(df[,2]+df[,3]+df[,4])/3,#28
									(df[,13]+df[,14])/2,#mor
									(df[,5]+df[,6])/2,#bla
									(df[,10]+df[,11]+df[,12])/3,#gas
									(df[,20]+df[,21])/2,#tro
									(df[,7]+df[,8]+df[,8])/3,#Dla
									(df[,17]+df[,18]+df[,19])/3	#spa
									)
colnames(methINTByStage)<-c("GeneID","ovo","x28","mor","bla","gas","tro","Dla","spa")


#jointure methINT et expression
int.exp<-merge (methINTByStage, mRNAByStage, by="GeneID")
int.exp<-as.data.frame(int.exp)
#filtrage expression nulle
int.exp<-int.exp[int.exp$expmoydev!=0,]


#Préparation du data.frame qui reçoit les résultats de corrélation
correlation<-data.frame(c(0),c(0),c(0),c(0))
colnames(correlation)<-c("id","p.value","r2","pente")
i<-1


#correlation
for (i in 1:nrow(int.exp))
{
#correlation[,1]<-exp.met
y<-as.numeric(t(int.exp[i,2:9]))
x<-as.numeric(t(int.exp[i,10:17]))
lm.x.y<-lm(y~x)
resultat<-summary(lm(y~x))[[4]]
p<-resultat[2,4]
r2<-summary(lm(y~x))[[8]]
pente<-resultat[2,1]
#correlation[i,1]<-id
#correlation[i,1]<-exp.met[i,1]
correlation[i,2]<-p
correlation[i,3]<-r2
correlation[i,4]<-pente
}
correlation[,1]<-int.exp[,1]

#write.table(correlation, "resultat.txt", sep=";", row.names=F,col.names=T)

#sélectionner les corrélations inférieures à un seuil
correlation.significative<-correlation[correlation$p.value<=0.05,]

#write.table(correlation.significative, "resultat_signif.txt", sep=";", row.names=F,col.names=T)

#représentation graphique
id.significatif<-row.names(correlation.significative)
id.significatif<-as.integer(id.significatif)

summary(correlation.significative)
summary(correlation.significative$pente>0)

plot(x=int.exp[id.significatif,]$ovo, y=int.exp[id.significatif,]$expov)
#save.image("evsgbMsovo.jpg")
plot(x=int.exp[id.significatif,]$x28, y=int.exp[id.significatif,]$exp2.8cells)
plot(x=int.exp[id.significatif,]$mor, y=int.exp[id.significatif,]$expmo)
plot(x=int.exp[id.significatif,]$bla, y=int.exp[id.significatif,]$expbl)
plot(x=int.exp[id.significatif,]$gas, y=int.exp[id.significatif,]$expga)
plot(x=int.exp[id.significatif,]$tro, y=int.exp[id.significatif,]$exptr)
plot(x=int.exp[id.significatif,]$Dla, y=int.exp[id.significatif,]$expdl)
plot(x=int.exp[id.significatif,]$spa, y=int.exp[id.significatif,]$expsp)


