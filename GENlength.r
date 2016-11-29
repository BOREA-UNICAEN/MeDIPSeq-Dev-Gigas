#SCRIPT GENE LENGTH
#Modification du répertoire de travail (comme R est un logiciel développé sous linux il faut utiliser des / au lieu des \)
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/Guillaume_SourceFiles") 
getwd()

GENlength<-read.table("Cgigas_v9_gene.gff", header=F, sep= "\t") #load feature gff file
GENlength<-data.frame(GENlength$V9, GENlength$V5-GENlength$V4) #compute feature length
colnames(GENlength)<-c("GeneID", "GENlengthgth") #add column names
head(GENlength)
GENlength$GeneID<-gsub("[ID=;]","",GENlength$GeneID) #remove crappy characters in GeneID column
GENlength$GeneID<-gsub("CG","CGI",GENlength$GeneID)
write.table(GENlength,"GENlength.txt", sep="\t", quote=F,col.names=T,row.names=F) #saves file