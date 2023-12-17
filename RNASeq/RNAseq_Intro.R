#RNA-seq Analysis Tutorial - Bioinformatics 2023

#Tutorial from https://combine-australia.github.io/RNAseq-R/ and
#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#Create an R Project to work in: https://swcarpentry.github.io/r-novice-gapminder/02-project-intro/index.html

#Download all data files from the following websites and put them in your project data folder
#Mouse mammary data (counts): https://figshare.com/s/1d788fd384d33e913a2a
#Drosophila data (counts): https://figshare.com/s/e08e71c42f118dbe8be6

#Install required packages (limma, edgeR, Glimma, org.Mm.eg.db, gplots, RcolorBrewer)
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

#Load each package independently to check for errors 
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

# Read the data into R - make sure to check file names
seqdata <- read.delim("data/GSE60450_LactationGenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)

head(seqdata)
dim(seqdata)
sampleinfo

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)

colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
head(countdata)

table(colnames(countdata)==sampleinfo$SampleName)

y <- DGEList(countdata)
# have a look at y
y

names(y)
y$samples

group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group

#change group to factor
group <- factor(group)
#take a look
group

# Add the group information into the DGEList
y$samples$group <- group
y$samples

columns(org.Mm.eg.db)
ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
# Have a look at the annotation
head(ann)
table(ann$ENTREZID==rownames(y$counts))

y$genes <- ann

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#CHALLENGE
# Second sample - Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5, h=10, col="blue")

y <- y[keep, keep.lib.sizes=FALSE]
y$samples$lib.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)

col.status <- c("blue","red","black")[sampleinfo$Status]
col.status

plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

#CHALLENGE
col.status <- c("lightblue", "pink", "lavender")[sampleinfo$Status]
pch.cell <- c(1, 4)[sampleinfo$CellType]

plotMDS(y, col=col.status, pch=pch.cell, cex=1.5) 
legend("topleft",fill=c("lightblue", "pink", "lavender"),legend=levels(sampleinfo$Status),cex=0.8)
legend("bottom",pch=c(1,4),legend=levels(sampleinfo$CellType),cex=0.8) 
title("Status & Cell Type")

# Dimension 3 appears to separate pregnant samples from the rest. Dim4?
plotMDS(y,dim=c(3,4),col=col.status,pch=pch.cell,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))

labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
glMDSPlot(y, labels=labels, groups=group, folder="mds")

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()

#CHALLENGE
#Heatmap Challenge 1
mypalette <- brewer.pal(11,"PiYG")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

#Heatmap Challenge 2
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row", labCol=group)

#Heatmap Challenge 3
low_select_var <- names(sort(var_genes))[1:500]
low_variable_lcpm <- logcounts[low_select_var,]
heatmap.2(low_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples

par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

par(mfrow=c(1,2))
plotMD(y,column = 7)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

#CHALLENGE
par(mfrow=c(1,2))
plotMD(logcounts,column = 11)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

sampleinfo
sampleinfo <- read.delim("data/SampleInfo_Corrected.txt", stringsAsFactors = TRUE)
sampleinfo
group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep="."))
y$samples$group <- group

# Look at group variable again
group

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design

## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

v
names(v)

y

#CHALLENGE
#1) It has the groups and normalization factors associated with each fly gene, this is the same as th y$samples data
#2) 15804 x 12

par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# Fit the linear model
#lmFit needs the voom object and the design matrix
fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate, levels=design)
cont.matrix

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#CHALLENGE - 2nd argument
cont.matrix2 <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate, L.PregVsLac = luminal.pregnant - luminal.lactate, levels=design)
cont.matrix2

fit.cont2 <- contrasts.fit(fit, cont.matrix2)
fit.cont2 <- eBayes(fit.cont2)
dim(fit.cont2)

summa.fit2 <- decideTests(fit.cont2)
summary(summa.fit2)

par(mfrow = c(1,1.7))
vennDiagram(summa.fit2,
            include= c("up", "down"),
            counts.col=c("purple", "black"),
            circle.col =c("purple", "green"))

# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.PregVsLac"], values = c(-1, 1), hl.col=c("blue","red"))

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="B.PregVsLac")

#CHALLENGE
volcanoplot(fit.cont2,coef=1,highlight=200,names=fit.cont$genes$SYMBOL, main="B.PregVsLac")


par(mfrow=c(1,3))
# Let's look at the first gene in the topTable, Wif1, which has a rowname 24117
stripchart(v$E["24117",]~group)
# This plot is ugly, let's make it better
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Wif1")

#topTable
topTable(fit.cont2)

#CHALLENGE
stripchart(v$E["12992",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Csn1s2b")

sampleinfo
sampleinfo <- read.delim("data/SampleInfo_Corrected.txt", stringsAsFactors = TRUE)
sampleinfo
group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep="."))

group2 <- group
levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate", "lum.preg", "lum.virgin")
glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", ylab="B", main="B.PregVsLac",
         counts=v$E, groups=group2, status=summa.fit[,1],
         anno=fit.cont$genes, side.main="ENTREZID", folder="volcano")

fit.treat <- treat(fit.cont2,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)

topTable(fit.treat,coef=1,sort.by="p")

# Notice that much fewer genes are highlighted in the MAplot
par(mfrow=c(1,2))
plotMD(fit.treat,coef=1,status=res.treat[,"B.PregVsLac"], values=c(-1,1), hl.col=c("blue","red"))
abline(h=0,col="grey")
plotMD(fit.treat,coef=2,status=res.treat[,"L.PregVsLac"], values=c(-1,1), hl.col=c("blue","red"))
abline(h=0,col="grey")

#CHALLENGE
fit.treat2 <- treat(fit.cont2,lfc=0.5)
res.treat2 <- decideTests(fit.treat2)
summary(res.treat2)

topTable(fit.treat2,coef=1,sort.by="p")

glMDPlot(fit.treat, coef=1, counts=v$E, groups=group2,
         status=res.treat, side.main="ENTREZID", main="B.PregVsLac",
         folder="md")

go <- goana(fit.cont2, coef="B.PregVsLac",species = "Mm")
topGO(go, n=10)

colnames(seqdata)

m <- match(rownames(fit.cont2),seqdata$EntrezGeneID)
gene_length <- seqdata$Length[m]
head(gene_length)

# Rerun goana with gene length information
go_length <- goana(fit.cont2,coef="B.PregVsLac",species="Mm",
                   trend=gene_length)
topGO(go_length, n=10)

#CHALLENGE
# Rerun goana with gene length information
go_length <- goana(fit.cont2,coef="L.PregVsLac",species="Mm",
                   trend=gene_length)
topGO(go_length, n=10)

# Load in the mouse c2 gene sets as Mm.c2
Mm.c2 <- load("~/Documents/School/Senior Year - 2022-23/bioinfo/RNASeq/data/mouse_c2_v5.rdata")
# Have a look at the first few gene sets
names(Mm.c2)[1:5]
# Number of gene sets in C2
length(Mm.c2)

c2.ind <- ids2indices(Mm.c2, rownames(v))
gst.camera <- camera(v,index=c2.ind,design=design,contrast = cont.matrix[,1],inter.gene.cor=0.05)
gst.camera[1:5,]
table(gst.camera$FDR < 0.05)

#CHALLENGE
gst.camera2 <- camera(v,index=c2.ind,design=design,contrast = cont.matrix[2,],inter.gene.cor=0.05)
gst.camera2[1:5,]
table(gst.camera2$FDR < 0.05)

Mm.H <- load("~/Documents/School/Senior Year - 2022-23/bioinfo/RNASeq/data/mouse_H_v5.rdata")
c2.ind2 <- ids2indices(Mm.H, rownames(v))
gst.camera3 <- camera(v,index=c2.ind2,design=design,contrast = cont.matrix[,1],inter.gene.cor=0.05)
gst.camera3[1:5,]
table(gst.camera3$FDR < 0.05)
