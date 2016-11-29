############## METHYLATION MAPS OF SCAFFOLDS #############
rm(list=ls(all=TRUE)) #empty memory
Sys.setenv(http_proxy="http://proxy.unicaen.fr:3128/") #proxy settings for R internet access

##########
setwd("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/methylmaps/bedfiles") 
getwd()



### Preparing data
## Chargement des données, originales dans Elizabeth Data\BAMtoBED  
message("loading data...")

#bed sample files
#df<-read.table("troc2s.bed", header=F, sep="\t")
#head(df)

# elimination des lignes ne correspondant pas au scaffold à mapper pour GenomeCoverage de BEDtools #
#df1<-df[df[,1]=="scaffold22",]
#head(df1)
#write.table(df1,"troc2s_scaffold22.bed", quote=F, col.names=F, row.names=F, sep="\t")

#The data *sample*_scaffold22.bed were processed using BEDtools GeneomeCoverage on Galaxy because my Laptop doesn't have a linux distribution installed...
#(bedtools_genomecoveragebed/2.22.0)	Output type 	hist 	
#Specify max depth 	0 	
#Genome file 	95: Scaffold22100bpintervals 	
#Treat split/spliced BAM or BED12 entries as distinct BED intervals when computing coverage. 	False 	
#Calculate coverage based on 	both strands combined 	
#Report the depth at each genome position with 1-based coordinates 	True 	
#Report the depth at each genome position with 0-based coordinatess 	False 	
#Calculate coverage of 5’ positions 	False 	
#Calculate coverage of 3’ positions 	False

# files were renamed *sample*_covsca22 and stored in the folder 'Buinages_Guillaume/methylmaps/Bedfiles/covsca22 files
#loading 1-base sample maps
#something like fileslist<-list.files("covsca22 files/") could have been smarter...
message("loading 1-base sample maps...")
ovo1<-read.table("covsca22 files/ovo1_covsca22.txt", header=F, sep="\t")
ovo2<-read.table("covsca22 files/ovo2_covsca22.txt", header=F, sep="\t")
ovo3<-read.table("covsca22 files/ovo3_covsca22.txt", header=F, sep="\t")
x281<-read.table("covsca22 files/28cell1_covsca22.txt", header=F, sep="\t")
x282<-read.table("covsca22 files/28cell2_covsca22.txt", header=F, sep="\t")
x283<-read.table("covsca22 files/28cell3_covsca22.txt", header=F, sep="\t")
mor1<-read.table("covsca22 files/morula1_covsca22.txt", header=F, sep="\t")
mor2<-read.table("covsca22 files/morula2_covsca22.txt", header=F, sep="\t")
bla2<-read.table("covsca22 files/blastula2_covsca22.txt", header=F, sep="\t")
bla3<-read.table("covsca22 files/blastula3_covsca22.txt", header=F, sep="\t")
gas1<-read.table("covsca22 files/gastrula1_covsca22.txt", header=F, sep="\t")
gas2<-read.table("covsca22 files/gastrula2_covsca22.txt", header=F, sep="\t")
gas3<-read.table("covsca22 files/gastrula3_covsca22.txt", header=F, sep="\t")
tro1<-read.table("covsca22 files/troc1_covsca22.txt", header=F, sep="\t")
tro2<-read.table("covsca22 files/troc2_covsca22.txt", header=F, sep="\t")
dla1<-read.table("covsca22 files/dlarv1_covsca22.txt", header=F, sep="\t")
dla2<-read.table("covsca22 files/dlarv2_covsca22.txt", header=F, sep="\t")
dla3<-read.table("covsca22 files/dlarv3_covsca22.txt", header=F, sep="\t")
spa1<-read.table("covsca22 files/spat1_covsca22.txt", header=F, sep="\t")
spa2<-read.table("covsca22 files/spat2_covsca22.txt", header=F, sep="\t")
spa3<-read.table("covsca22 files/spat3_covsca22.txt", header=F, sep="\t")

list<-list(ovo1=ovo1,ovo2=ovo2,ovo3=ovo3,x281=x281,x282=x282,x283=x283,mor1=mor1,mor2=mor2,bla2=bla2,bla3=bla3,gas1=gas1,gas2=gas2,gas3=gas3,tro1=tro1,tro2=tro2,dla1=dla1,dla2=dla2,dla3=dla3,spa1=spa1,spa2=spa2,spa3=spa3) # make a list with sample data tables

# loading nb reads per sample for normalization
read.number<-read.table("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/Guillaume_SourceFiles/readnumber.txt", header=T, sep="\t") # loads the nb_reads file
rownames(read.number)<-c("ovo1","ovo2","ovo3","x281","x282","x283","mor1","mor2","bla2","bla3","gas1","gas2","gas3","tro1","tro2","dla1","dla2","dla3","spa1","spa2","spa3") #gives a rowname corresponding to the sample name in the list

### formatting scaffold 22 physical methylation map per sample
for (i in 1:length(list))
{
colnames(list[[i]])<-c("chrom","base","cov") # gives column names to the objects in the list 
list[[i]]<-data.frame(list[[i]]$chrom,list[[i]]$base,((list[[i]]$cov/read.number[i,2])*1000)) # creates a data frame containing the base position and the normalised coverage per million read number
colnames(list[[i]])<-c("chrom","base","rpm") # gives normalized list objects the column names
list[[i]]<-list[[i]][1:(1000*(as.integer(nrow(list[[i]])/1000))),] # removes rows beyond the greatest 1000 multiple 
list[[i]]$levels<-dfcut<-cut(list[[i]]$base,as.integer(nrow(list[[i]])/10000)) # adds a 'levels' column and  cuts the base vector into intervals of 10kb
list[[i]]<-aggregate(list[[i]]$rpm, by=list(list[[i]]$levels),sum) # sums the read number (eg methylation) per interval 
colnames(list[[i]])<-c("levels","rpm") # gives normalized list objects the column names
list[[i]]$plotpos<-c(1:nrow(list[[i]])) # converts each interval into a position on the x axis ('plotpos') for plotting
}

### compute mean per stage Ovo, 28cells, I, spat and create data frame containing all the data for plotting
ovo<-data.frame(list[["ovo1"]]$plotpos, (list[["ovo1"]]$rpm )) #+ list[["ovo2"]]$rpm + list[["ovo3"]]$rpm)/3)
	colnames(ovo)<-c("plotpos","rpm.ovo")
x28<-data.frame(list[["x281"]]$plotpos, (list[["x281"]]$rpm + list[["x282"]]$rpm + list[["x283"]]$rpm)/3)
	colnames(x28)<-c("plotpos","rpm.x28")
I<-data.frame(list[["mor1"]]$plotpos, (	((list[["mor1"]]$rpm + list[["mor2"]]$rpm)/2 +
									 (list[["bla2"]]$rpm + list[["bla3"]]$rpm)/2 +
									 (list[["gas1"]]$rpm + list[["gas2"]]$rpm + list[["gas3"]]$rpm)/3 +
									 (list[["tro1"]]$rpm + list[["tro2"]]$rpm) /2 +
									 (list[["dla1"]]$rpm + list[["dla2"]]$rpm + list[["dla3"]]$rpm)/3)
									/5))
	colnames(I)<-c("plotpos","rpm.I")
spa<-data.frame(list[["spa1"]]$plotpos, (list[["spa1"]]$rpm + list[["spa2"]]$rpm + list[["spa3"]]$rpm)/3)
	colnames(spa)<-c("plotpos","rpm.spa")

stages<-list(ovo,x28,I,spa)
map<-Reduce(function(...) merge(..., all=T), stages) # builds a unique data table with the data
map[,2:5][map[,2:5]<=2]<-0 # Gives a 0 value for a sum of rpm <2 (considered background)

#write.table(map, "map.txt")


### PLOT
par(plot(map$rpm.spa~map$plotpos, 
			type="s",
			xlim=c(0,nrow(map)),
			ylim=c(0,max(map$rpm.spa)),
			col="blue",
			axes=F,
			ylab="methylation level",
			xlab=""),
			new=T)
x <- seq_along(map$rpm.spa)
y2 <- rep(map$rpm.spa, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(0,0,1,0.5)), new=T)
par(lines(x2, y2), new=T)

par(plot(map$rpm.I~map$plotpos, type="s",xlim=c(0,nrow(map)),ylim=c(0,max(map$rpm.spa)),col="orange", axes=F, xlab="",ylab=""), new=T)
x <- seq_along(map$rpm.I)
y2 <- rep(map$rpm.I, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(1,1,0,0.5)),new=T)
par(lines(x2, y2), new=T)

par(plot(map$rpm.x28~map$plotpos, type="s",xlim=c(0,nrow(map)),ylim=c(0,max(map$rpm.spa)),col="red",axes=F, xlab="",ylab=""), new=T)
x <- seq_along(map$rpm.x28)
y2 <- rep(map$rpm.x28, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(1,0,0,0.5)),new=T)
par(lines(x2, y2),new=T)

par(plot(map$rpm.ovo~map$plotpos, type="s",xlim=c(0,nrow(map)),ylim=c(0,max(map$rpm.spa)),col="black",axes=F, xlab="",ylab=""), new=T)
x <- seq_along(map$rpm.ovo)
y2 <- rep(map$rpm.ovo, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(0,1,0,0.5)), new=T)
par(lines(x2, y2),new=T)
par(segments(0,0,196.45),new=T)

### Gene map on scaffold 22
##load data
genemap<-read.table("D:/BOULOT/RECHERCHE/MeDIPseq/analyse_MeDIP/ElizabethCrowell/Bouinages_Guillaume/methylmaps/sca22_genemap.txt", header=T, sep="\t")
## create dataframe for plotting
gmap<-as.data.frame(seq(1:1964558)); colnames(gmap)<-c("base")
#builds 'list' containing vectors of the bases in scaffold 22 bases that are within a gene 
list<-list()  #creates an empty list
length(list)<-nrow(genemap) 
for (i in (1:nrow(genemap)))
{
list[[i]]<-c(genemap[i,2]:genemap[i,3]) #
}
ingene<-as.data.frame(unique(unlist(list))) # keeps unique values
gmap$ingene<-numeric(length(gmap$"base")) # creates a 'ingene' column in gmap giving a "0" value to each base
gmap$ingene[which(gmap$"base" %in% ingene[,1])]<-1 # attributes a "1" value for bases within genes


#plot
x<-gmap$base 
y<-gmap$ingene
df<-data.frame(x,y)
par(mfrow=c(2,1))
plot(y~x,	
	xlim=c(0,nrow(df)),
	ylim=c(0,max(df$y)), 
	type="s", 
	col="black", 
	axes=F, 
	xlab="", 
	ylab="Scaffold 22"
	)
segments(0,0,nrow(df),0)
segments(0,max(df$y),nrow(df),max(df$y))
segments(0,0,0,max(df$y))
segments(nrow(df),0,nrow(df),max(df$y))
x <- seq_along(df$y)
y2 <- rep(df$y, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col="black"), new=T)
par(lines(x2, y2),new=T)





############TESTS################
t<-data.frame(c("a","b","c","d"),c(1,2,3,4))
u<-data.frame(c("a","b","c","d"),c(11,22,33,44))
v<-data.frame(c("a","b","c","d"),c(111,222,333,444))
w<-data.frame(c("t","u","v"),c(1,10,100))
list<-list(t=t,u=u,v=v)
for (i in 1:length(list))
{
colnames(list[[i]])<-c("pos","cov")
list[[i]]<-data.frame(list[[i]]$pos,(list[[i]]$cov/w[i,2])*10)
colnames(list[[i]])<-c("pos","rpm")
}




x<-c(1:1042)
y<-rnorm(1042)
df<-data.frame(x,y)
colnames(df)<-c("x","y") #genere des données
df<-df[1:1000*(as.integer(nrow(df)/1000)),] #enlève les dernieres lignes au delà du plus grand multiple de 1000 
df$levels<-dfcut<-cut(df$x,as.integer(nrow(df)/100)) #ajoute une colonne 'levels' qui coupe le vecteur x en intervalles de 100
df$sum<-
df$mean<-ave(df$y,df$levels,FUN=mean) #computes mean by interval
df<-df[!duplicated(df$levels),3:4] #keeps only one line per interval and its mean
df$plotpos<-c(1:nrow(df)) # converts each interval into a position on the x axis

#plot
x<-df$plotpos
y<-df$mean
plot(x,y)
polygon(c(x,x[length(x)]), c(y, y[1]), col='red') 
ggplot(dt) + stat_density(aes(x=value,color=site),geom="line",position="dodge")+
  geom_ribbon(data=subset(gg,site=="site1" & x>q1),
              aes(x=x,ymax=y),ymin=0,fill="red", alpha=0.5)+
  geom_ribbon(data=subset(gg,site=="site2" & x<q2),
              aes(x=x,ymax=y),ymin=0,fill="blue", alpha=0.5)
			  

#gene map
#interval vector building
x<-c(3,7,12)
y<-c(4,10,13)
df<-data.frame(x,y)
z<-as.data.frame(seq(1:2000000))
list<-list()
length(list)<-nrow(df)
for (i in (1:nrow(df)))
{
list[[i]]<-c(df[i,1]:df[i,2])
}
ingene<-unique(unlist(list))
for (i in 1:nrow(z))
{
if (i %in% ingene) z$ingene[i]<-1
else z$ingene[i]<-0
}
			  
# gene map plot
x<-seq(1:10)
y<-c(0,0,1,0,1,1,1,0,1,0)
df<-data.frame(x,y)
plot(y~x,xlim=c(0,nrow(df)),ylim=c(0,max(df$y)), type="s", col="black", axes=F, xlab="", ylab="")
abline(h=c(0,1),v=c(0,nrow(df)))
x <- seq_along(df$y)
y2 <- rep(df$y, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(0.25,0.25,0.25,1)), new=T)
par(lines(x2, y2),new=T)

x<-c(3,4,7,10,12,13)
y<-c(1,0,1,0,1,0)
df<-data.frame(x,y)
plot(x,y,xlim=c(0,21),ylim=c(0,1), type="s", col="black", axes=F, xlab="", ylab="")
abline(h=c(0,1),v=c(0,nrow(20)))
x <- seq_along(df$y)
y2 <- rep(df$y, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)
par(polygon(x3, y3, border=NA, col=rgb(0.25,0.25,0.25,1)), new=T)
par(lines(x2, y2),new=T)