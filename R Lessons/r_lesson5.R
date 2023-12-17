#r lesson 5 
library(tidyverse)
library(readxl)

metadata <- read_excel(path="raw_data/baxter.metadata.xlsx",
                       col_types=c(sample = "text", fit_result = "numeric", Site = "text", Dx_Bin = "text",
                                   dx = "text", Hx_Prev = "logical", Hx_of_Polyps = "logical", Age = "numeric",
                                   Gender = "text", Smoke = "logical", Diabetic = "logical", Hx_Fam_CRC = "logical",
                                   Height = "numeric", Weight = "numeric", NSAID = "logical", Diabetes_Med = "logical",
                                   stage = "text")
)
metadata <- mutate(metadata, Height = na_if(Height, 0))
metadata <- mutate(metadata, Weight = na_if(Weight, 0))
metadata <- mutate(metadata, Site = recode(.x=Site, "U of Michigan"="U Michigan"))
metadata <- mutate(metadata, Dx_Bin = recode(.x=Dx_Bin, "Cancer."="Cancer"))
metadata <- mutate(metadata, Gender = recode(.x=Gender, "f"="female", "m"="male"))

metadata <- rename_all(.tbl=metadata, .funs=tolower)
metadata <- rename(.data=metadata,
                   previous_history=hx_prev,
                   history_of_polyps=hx_of_polyps,
                   family_history_of_crc=hx_fam_crc,
                   diagnosis_bin=dx_bin,
                   diagnosis=dx,
                   sex=gender)

metadata <- mutate(metadata, diagnosis = factor(diagnosis, levels=c("normal", "adenoma", "cancer")))

alpha <- read_tsv(file="raw_data/baxter.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, invsimpson, coverage)

meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

meta_alpha %>%
  + group_by(diagnosis) %>%
  + summarize(mean_shannon = mean(shannon), sd_shannon = sd(shannon), sum_shannon = sum(shannon), variance_shannon = var(shannon))

meta_alpha %>%
  + group_by(diagnosis, sex) %>%
  + summarize(mean_shannon = mean(shannon), sd_shannon = sd(shannon), sum_shannon = sum(shannon), variance_shannon = var(shannon))

meta_alpha %>%
  + group_by(sex) %>%
  + summarise(mean_age = mean(age), sd_age = sd(age))

meta_alpha %>%
  group_by(site) %>%
  summarize(mean_shannon = mean(shannon), sd_shannon = sd(shannon), N=n()) %>%
  arrange(desc(mean_shannon)) %>%
  head(n=1)

meta_alpha %>%
  group_by(site) %>%
  summarize(mean_shannon = mean(shannon), sd_shannon = sd(shannon), N=n()) %>%
  top_n(n=1, mean_shannon) %>%
  pull(site)

