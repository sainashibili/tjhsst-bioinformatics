library(tidyverse)
library(readxl)

# pcoa <- read_tsv(file="raw_data/baxter.braycurtis.pcoa.axes",
#                  col_types=cols(group=col_character())
# )
# 
# metadata <- read_excel(path="raw_data/baxter.metadata.xlsx",
#                        col_types=c(sample = "text", fit_result = "numeric", Site = "text", Dx_Bin = "text",
#                                    dx = "text", Hx_Prev = "logical", Hx_of_Polyps = "logical", Age = "numeric",
#                                    Gender = "text", Smoke = "logical", Diabetic = "logical", Hx_Fam_CRC = "logical",
#                                    Height = "numeric", Weight = "numeric", NSAID = "logical", Diabetes_Med = "logical",
#                                    stage = "text")
# )
# metadata <- mutate(metadata, Height = na_if(Height, 0))
# metadata <- mutate(metadata, Weight = na_if(Weight, 0))
# metadata <- mutate(metadata, Site = recode(.x=Site, "U of Michigan"="U Michigan"))
# metadata <- mutate(metadata, Dx_Bin = recode(.x=Dx_Bin, "Cancer."="Cancer"))
# metadata <- mutate(metadata, Gender = recode(.x=Gender, "f"="female", "m"="male"))
# 
# metadata <- rename_all(.tbl=metadata, .funs=tolower)
# metadata <- rename(.data=metadata,
#                    previous_history=hx_prev,
#                    history_of_polyps=hx_of_polyps,
#                    family_history_of_crc=hx_fam_crc,
#                    diagnosis_bin=dx_bin,
#                    diagnosis=dx,
#                    sex=gender)
# 
# metadata <- mutate(metadata, diagnosis = factor(diagnosis, levels=c("normal", "adenoma", "cancer")))
# 
# dir.create("processed_data", showWarnings=FALSE)
# write_tsv(x=metadata, path='processed_data/baxter.metadata.tsv')
# 
# metadata_pcoa <- inner_join(metadata, pcoa, by=c('sample'='group'))
# 
# ggplot(metadata_pcoa, aes(x=axis1, y=axis2, color=diagnosis)) +
#   geom_point(shape=19, size=2) +
#   scale_color_manual(name=NULL,
#                      values=c("black", "blue", "red"),
#                      breaks=c("normal", "adenoma", "cancer"),
#                      labels=c("Normal", "Adenoma", "Cancer")) +
#   coord_fixed() +
#   labs(title="PCoA of Bray-Curtis Distances Between Stool Samples",
#        x="PCo Axis 1",
#        y="PCo Axis 2") +
#   theme_classic()

alpha <- read_tsv(file="raw_data/baxter.groups.ave-std.summary",
                  col_types=cols(group=col_character())) %>%
                  filter(!method == 'std') %>%
                  select(group, sobs, shannon, invsimpson, coverage)
meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))
                