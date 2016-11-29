#SCRIPT CDS LENGTH
#Modification du répertoire de travail (comme R est un logiciel développé sous linux il faut utiliser des / au lieu des \)
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/Guillaume_SourceFiles") 
getwd()

CDSlength<-read.table("Cgigas_v9_exon.gff", header=F, sep= "\t") #load feature gff file
CDSlength<-data.frame(CDSlength$V9, CDSlength$V5-CDSlength$V4) #compute feature length
colnames(CDSlength)<-c("GeneID", "CDSlength") #add column names
head(CDSlength)
CDSlength$GeneID<-gsub("[Parent=;]","",CDSlength$GeneID) #remove crappy characters in GeneID column

#compute CDS length per gene
for (i in 2:nrow(CDSlength))
{
if (CDSlength$GeneID[i] == CDSlength$GeneID[i-1]) 
	CDSlength$CDSlength[i] <- CDSlength$CDSlength[i] + CDSlength$CDSlength[i-1]
else CDSlength$CDSlength[i] <- CDSlength$CDSlength[i] 
}#add CDSlength values if GeneID i =GeneID i-1

for (i in 2:nrow(CDSlength))
{
if (CDSlength$GeneID[i] == CDSlength$GeneID[i-1])
	CDSlength$CDSlength[i-1] <- 0
else CDSlength$CDSlength[i] <- CDSlength$CDSlength[i] 
}#replace double of i per 0 except the last

CDSlength<-CDSlength[CDSlength$CDSlength != 0,] #ne garde que le GeneID différent de 0
head(CDSlength)
write.table(CDSlength,"CDSlength.txt", sep="\t",quote=F, col.namnes=T, row.names=F) #saves file