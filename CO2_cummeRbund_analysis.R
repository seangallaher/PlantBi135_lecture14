# This is an R script to analyze 
# the RNA-Seq data from an analysis 
# of the carbon concentrating mechanism
# in Chlamydomonas reinhardtii as
# described in Fang et al., TPC 2012
# www.plantcell.org/cgi/doi/10.1105/tpc.112.097949

# For this analysis, we will be using
# a collection of functions 
# (i.e. a library) of tools specifically
# for analyzing RNA-Seq data. It is
# called cummeRbund, and we load it
# with the library() function. Notice 
# the capitol R in the middle. R is 
# case sensitive, so it matters that 
# all letters are lower case except for 
# that R. Here is the command: 

library(cummeRbund)

# Next, we need to load the RNA-Seq data.
# The initial analysis, normalizing the raw
# reads, has already been completed and 
# is ready to be loaded into a database
# for cummeRbund to work with. We will
# call this database "cuff", and we 
# will load the data into it with the
# readCufflinks() function like so:

cuff<-readCufflinks(dir = "RNAseq_data")

# To check that the data has been processed 
# correctly, we can type the name of the 
# database and it should return a brief
# summary of the data:

cuff

# You should see
# CuffSet instance with:
#	 6 samples 
#	 17741 genes ...


# Next, we will make a dendrogram of samples.
# This shows how closely related each sample
# is to the others. From this, we should
# be able to see if the CO2 concentration
# or the cia5 mutant has a larger effect
# on the overall transcriptome. Try it:

dend<-csDendro(genes(cuff))
dend

# Does the mutant or the CO2 concentration have
# a larger effect? 


# We may want to look at the expression
# of a specific gene. For example, HLA3 
# is important for the CCM in Chlamydomonas.
# The gene ID for HLA3 if Cre02.g097800.
# First, we need to get the expression
# data for that out of the cuff database:

myGeneData <- getGene(cuff, "Cre02.g097800")

# Next, we can see what its expression is
# with the expressionPlot() function:

expressionPlot(myGeneData)

# What do you see? When is HLA3 
# expressed the most?

# When we used expressionPlot, 
# we entered the gene ID for HLA3. 
# In R, it is preferable to use
# a variable for this so we
# can use the same code over and over again.
# Let's create a variable called
# myGeneID and put the HLA3 gene ID
# into it:

myGeneID<-"Cre02.g097800" 

# Now we can use expressionPlot()
# with our new variable:

expressionPlot(myGeneID)

# Now, we can change the value
# of myGeneID and reuse the same code.
# Try putting LCI1 (gene ID = Cre03.g162800)
# into myGeneID:

myGeneID<-"Cre03.g162800" 

# Go back and re-run the line that
# says expressionPlot(myGeneID)

# Try it with a few more genes:

myGeneID<-"Cre16.g662600" # LCI11C
myGeneID<-"Cre10.g436550" # LCI5 / EPYC1
myGeneID<-"Cre06.g284100" # RHP1


# A number of genes were known to code
# for proteins that are important for 
# the CCM. I have made a list of their
# gene IDs. We are going to create 
# a "vector" of these gene IDs with
# the c() function:

myGeneIds<-c("Cre02.g097800","Cre16.g662600",
	"Cre03.g162800","Cre04.g223300",
	"Cre06.g309000","Cre02.g112333",
	"Cre05.g236650","Cre06.g273950",
	"Cre04.g223100","Cre08.g367500",
	"Cre03.g151650","Cre09.g399552",
	"Cre07.g321750","Cre11.g477350",
	"Cre08.g367400","Cre04.g223250",
	"Cre10.g439700","Cre26.g756747",
	"Cre16.g651050","Cre16.g685000",
	"Cre05.g248450","Cre05.g248400",
	"Cre03.g204577","Cre12.g555700",
	"Cre07.g321800","Cre10.g426800",
	"Cre06.g281600")


# Make a mini database for just those genes:
myGenes<-getGenes(cuff,myGeneIds)

# Did it work? Check by looking for
# the summary info for myGenes:

myGenes

# cummeRbund has a built in fucntion
# for making heatmaps. Try it for your
# subset of genes:

csHeatmap(myGenes)

# Another useful function that cummeRbund
# has is called findSimilar(). You give it
# one gene ID, and it finds a number of genes
# with the most similar expression pattern.
# You can find as many similar genes as you
# want with the "n =" argument. Let's find
# the 5 genes whose expression patterns are
# most similar to HLA3:


mySimilar<-findSimilar(cuff,"Cre02.g097800",n = 5)

# Next, plot their expression:

mySimilar.expression<-expressionPlot(mySimilar,logMode=TRUE,showErrorbars=TRUE)

mySimilar.expression

# This data would probably be better as a heatmap.

mySimilar.heatmap<-csHeatmap(mySimilar,cluster='rows')
mySimilar.heatmap

##################################################

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



