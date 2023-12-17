#Saina Shibili
#making the lefse file

library(ANCOMBC)
library(BiocManager)
library(ComplexHeatmap)
library(DESeq2)
library(phyloseq)
library(tidyverse)
library(vegan)

#Making the phyloseq object
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

#sample.data.frame
sample.data.frame <- function(ps) {
  return(as(phyloseq::sample_data(ps), "data.frame"))
}

#otu.matrix
otu.matrix <- function(ps, force.taxa.cols = TRUE) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (taxa_are_rows(ps) & force.taxa.cols) {
    mat <- t(mat)
  }
  return(mat)
}

#taxa.matrix
taxa.matrix <- function(ps) {
  return(
    as(phyloseq::tax_table(ps), "matrix")
  )
}

#phyloseq2lefse
phyloseq2lefse <- function(
    ps,
    covars,
    file.name = "lefse_data.txt",
    taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    transpose.otus = TRUE
) 
  
  {
  if (length(covars) > 2) {
    warning(
      "The length of the `covars` vector is greater than 2. File must be edited manually for use in LEfSe, which throws and error when there are more rows than specific classes, subclasses, and subjects."
    )
  }
  smpl.data <- sample.data.frame(ps)
  smpl.data$Sample <- row.names(smpl.data)
  t.smpl.data <- t(smpl.data)
  t.smpl.data <- as.data.frame(t.smpl.data[c("Sample", covars), ])
  if (transpose.otus) {
    otu.tbl <- t(otu.matrix(ps)) # grab the otu table from the phyloseq object
  } else {
    otu.tbl <- otu.matrix(ps) # grab the otu table from the phyloseq object
  }
  tax.tbl <- taxa.matrix(ps) # grab the taxa table from the phyloseq object and coerce into a matrix
  tax.tbl <- tax.tbl[, taxa.levels]
  
  # The following loop goes through the taxa table starting from the highest level and gradually moving to the lowest taxonomic level. For each level it appends the unique taxonomic levels to the uniq.lvls vector.
  uniq.lvls <- c()
  for (i in c(1:length(tax.tbl[1, ]))) {
    lvls <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl[, 1:i]), sep = "|")))
    names(lvls) <- "tax.lvl"
    uniq.i <- as.character(unique(lvls$tax.lvl))
    uniq.lvls <- c(uniq.lvls, uniq.i)
  }
  tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep = "|")))
  row.names(tax.tbl.join) <- row.names(tax.tbl)
  names(tax.tbl.join) <- "tax.lvl"
  
  # This loop goes through each sample (which are now column names for t.smpl.data), and calculates the relative abundance for each unique taxonomic level (from above). These abundances only sum to 1 for *each taxonomic level*
  uniq.tax.lvl.abunds <- data.frame(row.names = uniq.lvls)
  for (smpl in names(t.smpl.data)) {
    abunds <- as.data.frame(otu.tbl[row.names(otu.tbl), smpl])
    total.abund <- sum(abunds[, 1])
    smpl.tax.lvl.abunds <- cbind(abunds, tax.tbl.join)
    
    smpl.uniq.lvl.abunds <- data.frame()
    for (uniq.lvl in uniq.lvls) {
      uniq.sub <- subset(smpl.tax.lvl.abunds, grepl(uniq.lvl, smpl.tax.lvl.abunds$tax.lvl, fixed = TRUE))
      lvl.total <- as.factor(sum(uniq.sub[, 1]) / total.abund)
      uniq.lvl.smpl <- data.frame(row.names = uniq.lvl, "sample" = lvl.total)
      names(uniq.lvl.smpl) <- smpl
      smpl.uniq.lvl.abunds <- rbind(smpl.uniq.lvl.abunds, uniq.lvl.smpl)
    }
    
    uniq.tax.lvl.abunds <- cbind(uniq.tax.lvl.abunds, smpl.uniq.lvl.abunds)
  }
  
  final.data <- rbind(t.smpl.data, uniq.tax.lvl.abunds)
  write.table(final.data, file=file.name, col.names=FALSE, sep="\t", quote=FALSE)
}

phyloseq2lefse(ps, "commons", "lefse_input.txt")
