#Saina Shibili
#r session 9

library(purrr)
library(broom)
library(tidyverse)

source('code/baxter.R')
alpha <- read_tsv(file="raw_data/baxter.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, invsimpson, coverage)
metadata <- get_metadata()
meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

taxonomy <- read_tsv(file="raw_data/baxter.cons.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")

otu_data <- read_tsv("raw_data/baxter.subsample.shared", col_types=cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(sample=Group) %>%
  pivot_longer(cols=-sample, names_to="otu", values_to="count") %>%
  mutate(rel_abund=count/10530)

source('code/baxter.R')
agg_phylum_data <- inner_join(otu_data, taxonomy) %>%
  group_by(sample, phylum) %>%
  summarize(agg_rel_abund=sum(rel_abund)) %>%
  inner_join(., get_metadata()) %>%
  ungroup() #without this, the sample and phylum columns remain grouped

top_phyla <- agg_phylum_data %>%
  group_by(phylum) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% # keep this so that the phyla are sorted properly
  top_n(5, median) %>%
  pull(phylum) # use pull to convert the names from a data frame to a vector of names

agg_phylum_data %>%
  filter(phylum %in% top_phyla) %>%
  mutate(phylum=factor(phylum, levels=top_phyla)) %>%
  ggplot(aes(x=phylum, y=agg_rel_abund, color=diagnosis)) +
  geom_boxplot() +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="There are no obvious phylum-level differences between the\ndiagnosis groups",
       x=NULL,
       y="Relative abundance") +
  theme_classic()

agg_phylum_data %>%
  filter(phylum %in% top_phyla) %>%
  mutate(phylum=factor(phylum, levels=top_phyla)) %>%
  ggplot(aes(x=phylum, y=agg_rel_abund, color=diagnosis)) +
  geom_jitter(pos=position_jitterdodge(jitter.width=0.2, dodge.width=0.8)) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="There are no obvious phylum-level differences between the\ndiagnosis groups",
       x=NULL,
       y="Relative abundance") +
  theme_classic()

source('code/baxter.R')

deep_taxonomy <- read_tsv(file="raw_data/baxter.cons.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="unclassified;", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=".*;", replacement=""))

agg_deep_data <- inner_join(otu_data, deep_taxonomy) %>%
  group_by(sample, taxonomy) %>%
  summarize(agg_rel_abund=sum(rel_abund)) %>%
  inner_join(., get_metadata()) %>%
  ungroup() #without this, the sample and phylum columns remain grouped

top_deep_taxa <- agg_deep_data %>%
  group_by(taxonomy) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% # keep this so that the phyla are sorted properly
  top_n(5, median) %>%
  pull(taxonomy) # use pull to convert the names from a data frame to a vector of names

agg_deep_data %>%
  filter(taxonomy %in% top_deep_taxa) %>%
  mutate(taxonomy=factor(taxonomy, levels=top_deep_taxa)) %>%
  mutate(agg_rel_abund=agg_rel_abund+1/21000) %>%
  ggplot(aes(x=taxonomy, y=agg_rel_abund, color=diagnosis)) +
  geom_boxplot() +
  geom_hline(yintercept=1/10530, color="gray") +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="There are no obvious phylum-level differences between the\ndiagnosis groups",
       x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  theme_classic()


