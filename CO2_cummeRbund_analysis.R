# This is an R script to analyze the RNA-Seq data from an analysis of the CCM in Fang et al., TPC 2012


# First we have to install the 
# cummeRbund library

BiocManager::install("cummeRbund")

# Load needed libraries
library(cummeRbund)

# Load the raw data into an object called "CO2data"
cuff<-readCufflinks(dir = "RNAseq_data")

# Check that the data has been processed correctly:
cuff

# Make a dendrogram of samples
dend<-csDendro(genes(cuff))
dend

# Make a dendrogram of samples with replicates
dend.rep<-csDendro(genes(cuff),replicates = TRUE)
dend.rep

# make a PCA plot of all
myPCA<-MDSplot(genes(cuff))
myPCA

# Make a volcano plot of all comparisons
#csVolcanoMatrix(genes(cuff))

# Check which samples are in the analysis
samples(genes(cuff))

# Make a volcano plot of specific comparisons
myVol<-csVolcano(genes(cuff), "WT_HiCO2", "cia5_VLCO2",alpha = 0.05,ylimits = c(0,2))
myVol

# Create a list of genes
myGeneIds<-c("Cre02.g097800","Cre16.g662600","Cre03.g162800","Cre04.g223300","Cre06.g309000","Cre02.g112333","Cre05.g236650","Cre06.g273950","Cre04.g223100","Cre08.g367500","Cre03.g151650","Cre09.g399552","Cre07.g321750","Cre11.g477350","Cre08.g367400","Cre04.g223250","Cre10.g439700","Cre26.g756747","Cre16.g651050","Cre16.g685000","Cre05.g248450","Cre05.g248400","Cre03.g204577","Cre12.g555700","Cre07.g321800","Cre10.g426800","Cre06.g281600")


# Make object for those genes
myGenes<-getGenes(cuff,myGeneIds)

# Make a heatmap for those genes:
csHeatmap(myGenes,cluster='both',replicates = TRUE)


# Plot one gene at a time:
myGeneID<-"Cre02.g097800" # HLA3
myGeneID<-"Cre16.g662600" # LCI11C
myGeneID<-"Cre03.g162800" # LCI1
myGeneID<-"Cre10.g436550" # LCI5 / EPYC1
myGeneID<-"Cre06.g284100" # RHP1

myGene<-getGene(cuff,myGeneID)

expressionPlot(myGene,logMode = FALSE)


# find significant genes
mySigGeneIds<-getSig(cuff,alpha=0.01,level='genes')

# How many significant genes are there?
length(mySigGeneIds)

# Create an object of only significant genes
mySigGenes<-getGenes(cuff,mySigGeneIds)




# Get FPKMs for significant genes
mySigFpkms<-fpkmMatrix(mySigGenes)


# Filter significant genes to remove any with less than 10 FPKMs
mySigMinGeneIds<-row.names(subset(mySigFpkms, (apply(mySigFpkms,1,max) > 10)))

# Create object with significant genes with >10 FPKMs
mySigMinGenes<-getGenes(cuff,mySigMinGeneIds)

# K-means clustering
myKclusters<-csCluster(mySigMinGenes, k = 16)

# plot K-means clusters
myKplot<-csClusterPlot(myKclusters)
myKplot

# Cluster 10 seems the most interesting. Isolate those gene IDs
myCluster<-subset(myKclusters$clustering, myKclusters$clustering == 10)

# How many are in this cluster
length(myCluster)

# What are their gene IDs
myClusterGeneIds<-names(myCluster)

# I get:
Cre01.g003950	Cre01.g003950
Cre02.g088551	Cre02.g088551
Cre02.g097800	HLA3
Cre03.g162800	LCI1
Cre03.g204577	DNJ31
Cre03.g212977	Cre03.g212977
Cre05.g237800	CYG64
Cre06.g281600	LCI23
Cre06.g309000	NAR1.2 (LCIA)
Cre07.g321750	Cre07.g321750
Cre07.g321800	Cre07.g321800
Cre08.g360550	ERM3
Cre08.g367500	LHCSR3.1
Cre09.g399552	LCR1
Cre10.g426800	Cre10.g426800
Cre11.g477350	Cre11.g477350
Cre12.g555700	DNAJ15
Cre16.g662600	LCI11C
Cre16.g685000	Cre16.g685000
Cre26.g756747	Cre26.g756747


# Find and plot 5 genes most similar to HLA3 Cre02.g097800


mySimilar<-findSimilar(cuff,"Cre02.g097800",n = 10)

mySimilar.expression<-expressionPlot(mySimilar,logMode=TRUE,showErrorbars=TRUE)

mySimilar.expression


# Plot as heatmap
mySimilar.heatmap<-csHeatmap(mySimilar,cluster='rows')
mySimilar.heatmap
