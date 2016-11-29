############## CpG density MAPS OF SCAFFOLDS #############
rm(list=ls(all=TRUE)) #empty memory
Sys.setenv(http_proxy="http://proxy.unicaen.fr:3128/") #proxy settings for R internet access

##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Elizabeth_Data/features_nbviewer") 
getwd()

### Preparing data
## Chargement des données, originales dans Elizabeth Data\BAMtoBED  
message("loading data...")

#bed sample files
df<-read.table("Cgigas_v9_CG.gff", header=F, sep="\t")
head(df)

#elimination des lignes ne correspondant pas au scaffold à mapper pour GenomeCoverage de BEDtools #
df1<-df[df[,1]=="scaffold22",]
df1<-data.frame(df1$"V1", df1$"V4", df1$"V5"); colnames(df1)<-c("chrom", "start","stop")
head(df1)
write.table(df1,"D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/methylmaps/CpG_scaffold22.txt", quote=F, col.names=F, row.names=F, sep="\t")

#Tool: Genome Coverage
#Name:	Genome Coverage on data 95 (scaffold 22 in 100 bp intervals) and data 163 (CpG_scaffold22.txt)
#Created:	Mon 29 Feb 2016 03:42:07 PM (UTC)
#Filesize:	38.3 MB
#Dbkey:	C.gigasv9
#Format:	tabular
#Galaxy Tool ID:	toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_genomecoveragebed/2.22.0
#Galaxy Tool Version:	2.22.0
#Tool Version:	bedtools v2.22.1
#Tool Standard Output:	stdout
#Tool Standard Error:	stderr
#Tool Exit Code:	0
#History Content API ID:	bbd44e69cb8906b559f2cc8aa501de5f
#Job API ID:	bbd44e69cb8906b5585ff04d328287ae
#History API ID:	1cc0139771d793b8
#UUID:	6f0497db-ab51-419f-a3e2-9162a66006de

#Input Parameter 	Value 	Note for rerun
#The BAM or BED file from which coverage should be computed 	163: CpG_scaffold22.txt 	
#Output type 	hist 	
#Specify max depth 	3 	
#Genome file 	95: Scaffold22100bpintervals 	
#Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage. 	False 	
#Calculate coverage based on 	both strands combined 	
#Report the depth at each genome position with 1-based coordinates 	True 	
#Report the depth at each genome position with 0-based coordinatess 	False 	
#Calculate coverage of 5’ positions 	False 	
#Calculate coverage of 3’ positions 	False
#######indicates whether the considered position bears a C or a G of a CpG dinucleotide

rm(list=ls(all=TRUE)) #empty memory
Sys.setenv(http_proxy="http://proxy.unicaen.fr:3128/") #proxy settings for R internet access

##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/methylmaps") 
getwd()

 

### formatting scaffold 22 CpG map: sums the number of CpG dinucleotides in windows of n bp (being successive x cordinates)giving a y value labeled as 'CpG density'
df<-read.table("CGmap_sca22[Galaxy170].tabular",header=F, sep="\t")
colnames(df)<-c("chrom", "pos", "CG") 
df$levels<-dfcut<-cut(df$pos,as.integer(nrow(df)/10000)) # adds a 'levels' column and cuts the 'pos'(ition on scaffold) vector into intervals of n bp (10000 here)
df<-aggregate(df$CG, by=list(df$levels),sum) # sums the number of CpGs per interval 
colnames(df)<-c("levels","CG") # gives column names or messes with same data.frame and variable names ;)
df$x<-c(1:nrow(df)) # converts each interval into a position on the x axis for plotting
df<-data.frame(df$"x",df$"CG"); colnames(df)<-c("x","y") # contains xy coordinates
plot(df$"y"~df$"x", type="l") 

