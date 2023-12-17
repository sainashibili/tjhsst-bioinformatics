#Saina Shibili
#R Lesson 7

#Activity 2
source("code/baxter.R")

pcoa <- read_tsv(file="raw_data/baxter.braycurtis.pcoa.axes",
                 col_types=cols(group=col_character())
                 )
metadata <- get_metadata()
metadata_pcoa <- inner_join(metadata, pcoa, by=c("sample"="group"))

ggplot(metadata_pcoa, aes(x=axis1, y=axis2, color=diagnosis)) +
  geom_point(shape=19, size=2) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  coord_fixed() +
  labs(title="PCoA of Bray-Curtis Distances Between Stool Samples",
       x="PCo Axis 1",
       y="PCo Axis 2") +
  theme_classic()

ggsave("code/figures/ordination.pdf")

rarefy <- read_tsv(file="raw_data/baxter.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-"))

#temperature of three diff cities on three diff days
temps <- tibble(day=c(1,2,3), chicago=c(75, 77, 74), detroit=c(69, 71, 70), nashville=c(79,80,78))
#tidys up temps
pivot_longer(temps, cols=c(-day), names_to='city', values_to='temperatures')

#rarefaction database
rarefy <- read_tsv(file="raw_data/baxter.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  pivot_longer(cols=c(-numsampled), names_to='sample', values_to='sobs') %>%
  mutate(sample=str_replace_all(sample, pattern="0.03-", replacement="")) %>%
  drop_na()

#joining rarefy and metadata
source('code/baxter.R')
metadata <- get_metadata()
metadata_rarefy <- inner_join(metadata, rarefy)

#roman metadata
roman_metadata <- metadata %>% 
                  mutate(stage=str_replace_all(stage, pattern="1", replacement="I")) %>%
                  mutate(stage=str_replace_all(stage, pattern="2", replacement="II")) %>%
                  mutate(stage=str_replace_all(stage, pattern="3", replacement="III")) %>%
                  mutate(stage=str_replace_all(stage, pattern="4", replacement= "IV"))
#counts how many of each stage there are
roman_metadata %>% count(stage)

#reduces number of samples
set.seed(1) #this makes sure that we all get the same result!
metadata_rarefy_sample <- metadata %>%
  sample_n(10) %>%
  inner_join(., rarefy)

#rarefaction curve
ggplot(metadata_rarefy_sample, aes(x=numsampled, y=sobs, group=sample, color=diagnosis, linetype=sex)) +
  geom_line() +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  scale_linetype_manual(name=NULL,
                        values=c(2, 6),
                        breaks=c("female", "male"),
                        labels=c("Female", "Male")) +
  coord_cartesian(xlim=c(0,20000), ylim=c(0,500)) +
  labs(title="Rarefaction curves are pretty pointless at this scale",
       x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()

#rarefaction curve with dots only at every 1000 samples
metadata_rarefy_sample_dots <- metadata_rarefy_sample %>% filter(numsampled %% 1000 == 0)

ggplot(metadata_rarefy_sample, aes(x=numsampled, y=sobs, group=sample, color=diagnosis)) +
  geom_line() +
  geom_point(data=metadata_rarefy_sample_dots, aes(x=numsampled, y=sobs, group=sample, color=diagnosis, shape=sex)) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  scale_shape_manual(name=NULL,
                     values=c(19, 1),
                     breaks=c("female", "male"),
                     labels=c("Female", "Male")) +
  coord_cartesian(xlim=c(0,20000), ylim=c(0,500)) +
  labs(title="Rarefaction curves are pretty pointless at this scale",
       x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()

#rarefaction curve with vertical line at 10530 samples
ggplot(metadata_rarefy_sample, aes(x=numsampled, y=sobs, group=sample, color=diagnosis)) +
  geom_line() +
  geom_vline(xintercept=10530) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  coord_cartesian(xlim=c(0,20000), ylim=c(0,500)) +
  labs(title="Rarefaction curves are pretty pointless at this scale",
       x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()

#scatter plot
source("code/baxter.R")
metadata <- get_metadata()

alpha <- read_tsv(file="raw_data/baxter.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, invsimpson, coverage)

meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

ggplot(meta_alpha, aes(x=bmi, y=shannon, color=diagnosis)) +
  geom_vline(xintercept=100, color="black") +
  geom_point(size=2) +
  scale_color_manual(name=NULL,
                     values=c("pink", "orchid1", "hotpink"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  labs(title="Shannon Diversity vs BMI",
       subtitle="Vertical line indicates the clinical screening threshold of 100",
       x="BMI",
       y="Shannon Diversity Index") +
  theme_classic()

#rarefaction curve
source("code/baxter.R")
metadata <- get_metadata()

rarefy <- read_tsv(file="raw_data/baxter.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  pivot_longer(cols=c(-numsampled), names_to="sample", values_to="sobs") %>%
  mutate(sample=str_replace_all(sample, pattern="0.03-", replacement="")) %>%
  drop_na()

metadata_rarefy <- inner_join(metadata, rarefy)

ggplot(metadata_rarefy, aes(x=numsampled, y=sobs, group=sample, color=diagnosis)) +
  geom_vline(xintercept=10530, color="gray", size=2) +
  geom_line() +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red"),
                     breaks=c("normal", "adenoma", "cancer"),
                     labels=c("Normal", "Adenoma", "Cancer")) +
  coord_cartesian(xlim=c(0,20000), ylim=c(0,500)) +
  labs(title="Rarefaction curve of patient samples",
       subtitle="Vertical line indicates the number of sequences that samples were rarefied to",
       x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()
ggsave("code/figures/rarefaction.pdf")


