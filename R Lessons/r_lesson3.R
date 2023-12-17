library(tidyverse)
library(readxl)

pcoa <- read_tsv(file="raw_data/baxter.braycurtis.pcoa.axes")
metadata <- read_excel(path="raw_data/baxter.metadata.xlsx",
                       col_types=c(sample = "text", fit_result = "numeric", Site = "text", Dx_Bin = "text",
                                   dx = "text", Hx_Prev = "logical", Hx_of_Polyps = "logical", Age = "numeric",
                                   Gender = "text", Smoke = "logical", Diabetic = "logical", Hx_Fam_CRC = "logical",
                                   Height = "numeric", Weight = "numeric", NSAID = "logical", Diabetes_Med = "logical",
                                   stage = "text")
)
metadata_pcoa <- inner_join(metadata, pcoa, by=c('sample'='group'))

ggplot(metadata_pcoa, aes(x=axis1, y=axis2, color=dx)) +
  geom_point(shape=19, size=2) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red", "black"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  coord_fixed() +
  labs(title="PCoA of Bray-Curtis Distances Between Stool Samples",
       x="PCo Axis 1",
       y="PCo Axis 2") +
  theme_classic()

ggsave("ordination.pdf")
