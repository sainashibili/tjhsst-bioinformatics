#r lesson 6

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

#calculating mean and sd of diagnoses
meta_alpha %>%
  group_by(diagnosis) %>%
  summarize(mean=mean(fit_result), sd=sd(fit_result))

#strip chart of fit and diagnosis
ggplot(meta_alpha, aes(x=diagnosis, y=fit_result, color=diagnosis)) +
  geom_jitter(shape=19, size=2) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  scale_x_discrete(limits=c("normal", "adenoma", "cancer"),
                   labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="Relationship between FIT result and subject's diagnosis",
       x=NULL,
       y="FIT Result") +
  theme_classic()

#strip chart that shows shannon diversity for each diagnosis category
ggplot(meta_alpha, aes(x=diagnosis, y=shannon, color=diagnosis)) +
    geom_jitter(shape=19, size=2) +
    scale_color_manual(name=NULL,
                       values=c("black", "blue", "red"),
                       breaks=c("normal", "adenoma", "cancer"),
                       labels=c("Normal", "Adenoma", "Cancer")) +
    scale_x_discrete(limits=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
    labs(title="Relationship between Shannon Diversity and subkect's diagnosis",
         x=NULL,
         y="Shannon Diversity") +
    theme_classic() + 
    coord_cartesian(ylim=c(0,5))

#box plot of fit and subject diagnosis
ggplot(meta_alpha, aes(x=diagnosis, y=fit_result, color=diagnosis)) +
  geom_boxplot() +
  #geom_boxplot(notch=TRUE) will give a notched boxplot
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  scale_x_discrete(limits=c("normal", "adenoma", "cancer"),
                   labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="Relationship between FIT result and subject's diagnosis",
       x=NULL,
       y="FIT Result") +
  theme_classic()

#violin plot of fit and subject diagnosis
ggplot(meta_alpha, aes(x=diagnosis, y=fit_result, fill=diagnosis)) +
  geom_violin() +
  scale_fill_manual(name=NULL,
                    values=c("black", "blue", "red"),
                    breaks=c("normal", "adenoma", "cancer"),
                    labels=c("Normal", "Adenoma", "Cancer")) +
  scale_x_discrete(limits=c("normal", "adenoma", "cancer"),
                   labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="Relationship between FIT result and subject's diagnosis",
       x=NULL,
       y="FIT Result") +
  theme_classic()

ggplot(meta_alpha, aes(x=diagnosis, y=shannon)) +
  geom_density_ridges()
