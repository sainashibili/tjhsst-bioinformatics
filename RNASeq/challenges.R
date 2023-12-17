#Saina Shibili
#RNASeq Challenges

#CHALLENGE 1
# Second sample - Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5, h=10, col="blue")

#CHALLENGE 2
col.status <- c("lightblue", "pink", "lavender")[sampleinfo$Status]
pch.cell <- c(1, 4)[sampleinfo$CellType]

plotMDS(y, col=col.status, pch=pch.cell, cex=1.5) 
legend("topleft",fill=c("lightblue", "pink", "lavender"),legend=levels(sampleinfo$Status),cex=0.8)
legend("bottom",pch=c(1,4),legend=levels(sampleinfo$CellType),cex=0.8) 
title("Status & Cell Type")

#CHALLENGE 3
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

#CHALLENGE 4
#1) It has the groups and normalization factors associated with each fly gene, this is the same as th y$samples data
#2) 15804 x 12

#CHALLENGE 5
stripchart(v$E["12992",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Csn1s2b")

#CHALLENGE 6
fit.treat2 <- treat(fit.cont2,lfc=0.5)
res.treat2 <- decideTests(fit.treat2)
summary(res.treat2)

topTable(fit.treat2,coef=1,sort.by="p")

glMDPlot(fit.treat, coef=1, counts=v$E, groups=group2,
         status=res.treat, side.main="ENTREZID", main="B.PregVsLac",
         folder="md")

#CHALLENGE 7
# Rerun goana with gene length information
go_length <- goana(fit.cont2,coef="L.PregVsLac",species="Mm",
                   trend=gene_length)
topGO(go_length, n=10)

#CHALLENGE
gst.camera2 <- camera(v,index=c2.ind,design=design,contrast = cont.matrix[2,],inter.gene.cor=0.05)
gst.camera2[1:5,]
table(gst.camera2$FDR < 0.05)

#CHALLENGE 8
Mm.H <- load("~/Documents/School/Senior Year - 2022-23/bioinfo/RNASeq/data/mouse_H_v5.rdata")
c2.ind2 <- ids2indices(Mm.H, rownames(v))
gst.camera3 <- camera(v,index=c2.ind2,design=design,contrast = cont.matrix[,1],inter.gene.cor=0.05)
gst.camera3[1:5,]
table(gst.camera3$FDR < 0.05)
