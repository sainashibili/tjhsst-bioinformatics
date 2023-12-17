#R lesson 2

#reading libraries
library(tidyverse)
library(readxl)

#reading and joining files
pcoa <- read_tsv(file="raw_data/baxter.braycurtis.pcoa.axes")
metadata <- read_excel(path="raw_data/baxter.metadata.xlsx")
metadata_pcoa <- inner_join(metadata, pcoa, by=c('sample'='group'))

#plotting data from files
ggplot(metadata_pcoa, aes(x=axis1, y=axis2, color=dx, shape=dx)) +
  geom_point(size=2) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  scale_shape_manual(name=NULL,
                     values=c(15, 16, 17),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  coord_fixed() +
  labs(title="PCoA of Bray-Curtis Distances Between Stool Samples",
       x="PCo Axis 1",
       y="PCo Axis 2") +
  theme_classic()

#saves as a pdf
ggsave("ordination.pdf")