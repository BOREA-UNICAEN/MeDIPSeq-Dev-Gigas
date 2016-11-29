############## mRNA expression VS CDS/INT METHYLATION RATIO DISTANCE FROM 1  #############
##### This scripts computes the length of genebody features (CDS and INT) per gene, then normalizes the 
# methylation level in INT vs CDS for each gene regarding INT and CDS length, respectively. 
# Then it computes the normalized INT/CDS methylation ratio per sample, then the mean by stage.
# Performs useless anova because we keep all the genes regardless meth pattern changes across development
# (but the user is free to play with that, although the overall outcome is always the same).
# Loads expression levels, merges tables, keeps column with mean mRNA levels. Excludes expression levels outlines
# (ie not expressed or with mean >20000 rpkm (n=8) which screws everything up). Merges both.Plots mRNA levels
# against ratio meth INT/CDS quantiles. Tests differences in mRNA level per quantile using ANOVA, pairwise student's t-test and linear model correlation (the user needs to deal with Inf values first for that last test).   


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
#feature length for normalisation for GFF files (oyster V9)
int <- read.table("Cgigas_v9_intron.gff", header = FALSE, sep = "")
cds <- read.table("Cgigas_v9_exon.gff", header = FALSE, sep = "")
gen <- read.table("Cgigas_v9_gene.gff", header = FALSE, sep = "")
#methylation
methCDS<-read.table("countsCDSfiltered.txt", header=T, sep="")
methINT<-read.table("countsINTfiltered.txt", header=T, sep="")
methGB<-read.table("countsGENfiltered.txt", header=T, sep="")
FractInCDS<-read.table("FractionInCDSDevAll.txt", header=T, sep="\t")
methCDSvsINT <- read.table("countsnCDS_INT_filtered_devall.txt",header = TRUE, sep="\t")
#Gene subsets
DevGiga<-read.table("gigaton_dev.txt",header = TRUE, sep="\t")
hox<-read.table("Homeobox_genes.txt",header = TRUE, sep="\t")

###### Compute normalized methINT/methCDS ratio 
## compute feature length by gene
# format INT len data data table
INTlen <- data.frame (int[,9], int[,5]-int[,4])
colnames (INTlen) <- c ("GeneID", "length")
INTlen$GeneID <-  substr(as.character(INTlen$GeneID), 8, nchar(as.character(INTlen$GeneID))-1)
# sums INT length for each gene
INTlen <- aggregate(INTlen$length, by = list(INTlen$GeneID), FUN=sum)
colnames (INTlen) <- c ("GeneID", "INTlength")

#does the same for CDS
CDSlen <- data.frame (cds[,9], cds[,5]-cds[,4])
colnames (CDSlen) <- c ("GeneID", "length")
CDSlen$GeneID <-  substr(as.character(CDSlen$GeneID), 8, nchar(as.character(CDSlen$GeneID))-1)
CDSlen <- aggregate(CDSlen$length, by = list(CDSlen$GeneID), FUN=sum)
colnames (CDSlen) <- c ("GeneID", "CDSlength")

#does the same for GEN
GENlen <- data.frame (gen[,9], gen[,5]-gen[,4])
colnames (GENlen) <- c ("GeneID", "length")
GENlen$GeneID <-  substr(as.character(GENlen$GeneID), 8, nchar(as.character(GENlen$GeneID))-1)
GENlen <- aggregate(GENlen$length, by = list(GENlen$GeneID), FUN=sum)
colnames (GENlen) <- c ("GeneID", "GENlength")

## builds table with teh results
FeatLen <- merge(CDSlen,INTlen, by = "GeneID", all = TRUE)
dim(FeatLen)


### Normalise methylation counts by feature length
mCDSn <- data.frame(methCDS$GeneID, methCDS[,-1]/CDSlen$CDSlength)
colnames (mCDSn) <- colnames (methCDS)
mINTn <- data.frame(methINT$GeneID, methINT[,-1]/INTlen$INTlength)
colnames (mINTn) <- colnames (methINT)
# combine in one data frame
methCDSvsINT <- merge (mCDSn, mINTn, by = "GeneID", all=TRUE)


#replaces meth INT NA by 0 meth CDS = 0 by 10-15 to have the script work (avoids /o)
# at the cost of introducing a 10-13% error in measuring the pattern, considered negligible.
#
methCDSvsINT[,2:41][is.na(methCDSvsINT[,2:41]) ] <- 0
methCDSvsINT[,2:21][methCDSvsINT[,2:21] == 0] <- 0.0000000000000001

# compute ratio INT/CDS meth per gène, orders by stage
ratioINT_CDS_DevAll <- data.frame(methCDSvsINT[,1],
							(methCDSvsINT[,35]/methCDSvsINT[,15]),
							(methCDSvsINT[,36]/methCDSvsINT[,16]), #ovo
							(methCDSvsINT[,22]/methCDSvsINT[,2]),
							(methCDSvsINT[,23]/methCDSvsINT[,3]), 
							(methCDSvsINT[,24]/methCDSvsINT[,4]), #28C 
							(methCDSvsINT[,33]/methCDSvsINT[,13]),
							(methCDSvsINT[,34]/methCDSvsINT[,14]),	#mor
							(methCDSvsINT[,25]/methCDSvsINT[,5]),
							(methCDSvsINT[,26]/methCDSvsINT[,6]), #bla
							(methCDSvsINT[,30]/methCDSvsINT[,10]),
							(methCDSvsINT[,31]/methCDSvsINT[,11]),
							(methCDSvsINT[,32]/methCDSvsINT[,12]), #gas
							(methCDSvsINT[,40]/methCDSvsINT[,20]),
							(methCDSvsINT[,41]/methCDSvsINT[,21]), #troc
							(methCDSvsINT[,27]/methCDSvsINT[,7]),
							(methCDSvsINT[,28]/methCDSvsINT[,8]), 
							(methCDSvsINT[,29]/methCDSvsINT[,9]), #Dla
							(methCDSvsINT[,37]/methCDSvsINT[,17]),
							(methCDSvsINT[,38]/methCDSvsINT[,18]),
							(methCDSvsINT[,39]/methCDSvsINT[,19]) #spa
						)
colnames(ratioINT_CDS_DevAll)<-c("GeneID","ovo1","ovo2", "x281" ,"x282","x283","mor1", "mor2","bla2","bla3", "gas1","gas2", "gas3","tro1","tro2","Dla1","Dla2","Dla3","spa1","spa2" ,"spa3")
head (ratioINT_CDS_DevAll)	


######### ANOVA ratioINTCDS vs DEV STAGE ###########
#ANOVA by row (factor is "development stage"):
ds <- ratioINT_CDS_DevAll
rownames(ds) <- ds[ ,1]; ds <- ds[ ,-1]
colnames(ds) <- substr(colnames(ds), 1, nchar(colnames(ds))-1)
testRes <- numeric(nrow(ratioINT_CDS_DevAll))

message("Running ANOVA against development stages...")

for (i in 1:nrow(ds))
{
 values <- as.numeric(ds[i, ])
 test <- aov(values ~ as.factor(colnames(ds)))
 testRes[i] <- summary(test)[[1]]$"Pr(>F)"[1]
}
ratioINT_CDS_DevAll<-data.frame(ratioINT_CDS_DevAll, testRes)

### Keep all genes regardless of ANOVA results 
df<-ratioINT_CDS_DevAll

# Exclude all genes rows for which ANOVA gave non-significant results:
#df<-ratioINT_CDS_DevAll[!is.na(ratioINT_CDS_DevAll$testRes) & ratioINT_CDS_DevAll$testRes< 0.01,] 
message(paste(dim(df)[1], " genes kept", sep = ""))

#Gene subsets
#df<-merge(hox, df,by="GeneID"); colnames(df)<-c(colnames(df))


#moyenne par stade
ovo<-numeric(nrow(df))
x28<-numeric(nrow(df))
mor<-numeric(nrow(df))
bla<-numeric(nrow(df))
gas<-numeric(nrow(df))
tro<-numeric(nrow(df))
Dla<-numeric(nrow(df))
spa<-numeric(nrow(df))
for (i in 1:nrow(df))
{
o<-c(df[i,2],df[i,3])
ovo[i]<-mean(o,na.rm=TRUE)
x<-c(df[i,4],df[i,5], df[i,6])
x28[i]<-mean(x,na.rm=TRUE)
m<-c(df[i,7],df[i,8])
mor[i]<-mean(m,na.rm=TRUE)
b<-c(df[i,9],df[i,10])
bla[i]<-mean(b,na.rm=TRUE)
g<-c(df[i,11],df[i,12],df[i,13])
gas[i]<-mean(g,na.rm=TRUE)
tr<-c(df[i,14],df[i,15])
tro[i]<-mean(tr,na.rm=TRUE)
d<-c(df[i,16],df[i,17],df[i,18])
Dla[i]<-mean(d,na.rm=TRUE)
s<-c(df[i,19],df[i,20],df[i,21])
spa[i]<-mean(s,na.rm=TRUE)
}

ratioINT_CDS_ByStage<-data.frame(df[,1],ovo,x28,mor,bla,gas,tro,Dla,spa,df$testRes)
colnames(ratioINT_CDS_ByStage)<-c("GeneID","ovo","x28","mor","bla","gas","tro","Dla","spa","testRes")
head(ratioINT_CDS_ByStage)


#merge ration and mRNA expresion per genes, keep only methylated genes
df1<-merge(ratioINT_CDS_ByStage,mRNAByStage,by="GeneID", all=FALSE)
#write.table(df1, "ratio_vs_mRNA.txt")
#exclusion expression aberrante 
#df1<-df1[df1$expmoydev!=0,]
df1<-df1[df1$expmoydev<20000,]
# computes 'distance from 1' for the ratio variable: represents skew from balanced methylation between CDS and INT
#df1[,2:9][df1[,2:9]==0] <- 0.01
df1[,2:9][df1[,2:9]<1] <- 1/(df1[,2:9][df1[,2:9]<1])

## moyennes ratio sur ensemble du dev 
df2<-data.frame(df1[1], 
					(df1[,2]+df1[,3]+df1[,4]+df1[,5]+df1[,6]+df1[,7]+df1[,8]+df1[,9])/8,
					df1[19]
				)
colnames(df2)<-c("GeneID", "ratiomean", "expmoydev")



#exclusion ratio=0 or Inf for correlation test 
#df2<-df2[df2$ratiomean!=0,] 
#df2<-df2[df2$ratiomean!=Inf,]


# representation graphique
#division en intervalles, breaks = nb d'intervalles désirés +1
df2$breaks <- cut( 	df2$ratiomean,
							#breaks= c(1,5,50,Inf),
							quantile(	df2$ratiomean, 
										na.rm = TRUE,
										prob = seq(0, 1, length = 4),
										type = 5),
							include.lowest = TRUE)

# visualisation des 6 premières lignes 
head(df2)

#définitinn de la variable quantiles et visualisation des valeurs
breaks<-as.factor(df2$breaks)
summary(breaks)													
									
						
													
#Représentation graphique en Boxplot de mRNACV ~ breaks methylation balance
boxplot(df2$expmoydev~breaks,
		range = 0.5, 
		width = NULL, 
		varwidth = F, # surface des plots en fonction de la racine carrée de n
		notch = T, # dessine les notch des boxes (svt signif si pas overlap)
		outline = FALSE, # ne trace pas les points sur le graphe
		col = c("red","blue","blue"),#rainbow(8), 
		log = "", # echelle log des axes log = "x,y"
		# le vecteur names doit faire la meme longieir que le nb de quantiles 
		names = c("1st","2nd","3rd"), 
		las = 1,
		xlab = "ratio INT/CDS (quantiles)",
		ylab = "mRNA level",
		main = "", 
		#ylim = c(-0.1,20)
		)

 summary(aov(df2$expmoydev~breaks))
 pairwise.t.test(df2$expmoydev, breaks, p.adj = "none")
 summary(lm(df2$expmoydev~df2$ratiomean))
 
