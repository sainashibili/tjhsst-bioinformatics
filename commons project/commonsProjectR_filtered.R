#Saina Shibili
#Commons project - filtered data

library(ANCOMBC)
library(BiocManager)
library(ComplexHeatmap)
library(DESeq2)
library(phyloseq)
library(tidyverse)
library(vegan)

otu <- read.table(file = "filtered-feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
taxonomy <- read.table(file = "filtered-taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Greengenes format
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

#cleaning up the tax file
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("filtered-tree.nwk")

ps <- phyloseq(OTU, TAX, SAMPLE,TREE)

#rarefaction
set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=2000, replace=F)

#COMMONS ANALYSIS
#alpha diversity plots for commons
plot_richness(ps.rarefied, x="commons", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

#pairwise test with Shannon rank-sum test for commons
rich = estimate_richness(ps.rarefied, measures = c("Observed", "Shannon"))
wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ps.rarefied)$commons, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.shannon

#beta diversity plots for commons - bray curtis/PCoA
dist = phyloseq::distance(ps.rarefied, method="bray")
ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)
plot_ordination(ps.rarefied, ordination, color="commons") + 
  theme_classic() +
  theme(strip.background = element_blank())

#permanova/adonis beta diversity for commons
metadata_adonis <- data.frame(sample_data(ps.rarefied))
test.adonis <- adonis2(dist ~ commons, data = metadata_adonis)
test.adonis

#binary jaccard beta diversity for commons - NMDS
dist = phyloseq::distance(ps.rarefied, method="jaccard", binary = TRUE)
ordination = ordinate(ps.rarefied, method="NMDS", distance=dist)
plot_ordination(ps.rarefied, ordination, color="commons") + 
  theme_classic() +
  theme(strip.background = element_blank())

#ANOSIM for commons 
metadata_anosim <- data.frame(sample_data(ps.rarefied))
anosim(dist, metadata_anosim$commons)

#ABUNDANCE BAR PLOTS

#phylum
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

#stricly by phylum -- not sorted by groups
ggplot(ps.melt, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))

ps.melt <- ps.melt %>%
  group_by(commons, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,commons,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum), width=0.8) + 
  labs(x="", y="%") +
  facet_wrap(~commons, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))

#genus
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Genus <- as.character(ps.melt$Genus)

#strictly by genus, not sorted by groups
ggplot(ps.melt, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus), width = 0.05) + 
  labs(x="", y="%") +
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))

ps.melt <- ps.melt %>%
  group_by(commons, Genus) %>%
  mutate(median=median(Abundance))
# select group mean > 1
keep <- unique(ps.melt$Genus[ps.melt$median > 2.5])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 2.5%"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,commons,Genus) %>%
  summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus)) + 
  labs(x="", y="%") +
  facet_wrap(~commons, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))

#DIFFERENTIAL ABUNDANCE - DESeq2
sample_data(ps)$commons <- as.factor(sample_data(ps)$commons) # factorize for DESeq2
ps.taxa.sub <- tax_glom(ps, taxrank = 'Species', NArm = FALSE)
#filter sparse features, with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)
ps_ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ commons)
# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)

matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Species"])
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
# Define the annotation color for columns and rows
annotation_col = data.frame(
#  Subject = as.factor(metadata_sub$subject), 
  `Commons` = as.factor(metadata_sub$commons), 
  check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)
rownames(annotation_row) = rownames(matrix)

#ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
#  Subject = c(`subject-1` = "red", `subject-2` = "blue"),
  `Commons` = c(Gandhi = "purple", Turing = "yellow", Hopper = "red", Carver = "blue", Newton = "green", Galileo = "orange"),
  Phylum = phylum_col
)

ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors)

#
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")
