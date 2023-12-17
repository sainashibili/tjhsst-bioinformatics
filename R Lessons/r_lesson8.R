#Saina Shibili
#r lesson 8 

#installing libraries
library(purrr)
library(broom)

#meta_alpha database
source('code/baxter.R')
alpha <- read_tsv(file="raw_data/baxter.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, invsimpson, coverage)
metadata <- get_metadata()
meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

#makes the summary_data table
meta_alpha %>%
  nest(data = -diagnosis) %>%
  mutate(summary_data=map(data, ~summary(.x$shannon) %>% tidy), 
         N = map(data, ~nrow(.x))) %>%
  unnest(cols= c(summary_data, N)) %>%
  select(-data)

#normal quantile plot
ggplot(meta_alpha, aes(sample=shannon, group=diagnosis, color=diagnosis)) + geom_qq() + stat_qq_line()

#transforming the data
meta_alpha <- mutate(meta_alpha, scaled_shannon=shannon^3)

#new normal quantile plot
ggplot(meta_alpha, aes(sample=scaled_shannon, group=diagnosis, color=diagnosis)) +
  geom_qq() + stat_qq_line()

#histograms
ggplot(meta_alpha, aes(x=shannon)) + geom_histogram()
ggplot(meta_alpha, aes(x=scaled_shannon)) + geom_histogram()

#checks normality -- old data
meta_alpha %>% pull(shannon) %>% shapiro.test()
#checks normality -- new data
meta_alpha %>% pull(scaled_shannon) %>% shapiro.test()

#aov test
diagnosis_shannon_aov <- aov(scaled_shannon~diagnosis, data=meta_alpha)
summary(diagnosis_shannon_aov)

#EXAMPLE - if experimental p-value less than 0.05 then use this, p-value is not less than 0.05
TukeyHSD(diagnosis_shannon_aov)

#test we would use if we hadnt transformed the values
result <- kruskal.test(shannon~diagnosis, data=meta_alpha)

#shows computer version of methods
glimpse(result)

#extracting the type of test
result[["test"]]
result$test

#determines which groups have a statistically significant difference
pairwise.wilcox.test(g=meta_alpha$diagnosis, x=meta_alpha$shannon, p.adjust.method="BH")

#significant difference in the number of OTUS by diagnosis group
ggplot(meta_alpha, aes(sample=sobs, group=diagnosis, color=diagnosis)) + geom_qq() + stat_qq_line()
meta_alpha <- mutate(meta_alpha, scaled_sobs=sobs^0.5)
ggplot(meta_alpha, aes(sample=scaled_sobs, group=diagnosis, color=diagnosis)) +
  geom_qq() + stat_qq_line()
ggplot(meta_alpha, aes(x=sobs)) + geom_histogram()
ggplot(meta_alpha, aes(x=scaled_sobs)) + geom_histogram()
diagnosis_sobs_aov <- aov(scaled_sobs~diagnosis, data=meta_alpha)
summary(diagnosis_sobs_aov)

#testing multiple hypothesis
meta_alpha %>%
  mutate(diagnosis = as.character(diagnosis),
         sex = as.character(sex), #unnecessary since it's already a character vector
         smoke = as.character(smoke)) %>%
  select(sample, shannon, diagnosis, sex, smoke) %>%
  pivot_longer(cols=c(diagnosis, sex, smoke), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(shannon ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

#sees whether variation in fit_result data is significant across diagnosis groups
meta_alpha %>%
  select(sample, fit_result, diagnosis, site) %>%
  nest(data = -site) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(fit_result ~ diagnosis, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

#correlation between fit_result and shannon
cor.test(meta_alpha$fit_result, meta_alpha$shannon)

#regression line
lm_shannon_bmi <- lm(shannon~bmi + diagnosis, data=meta_alpha)
summary(lm_shannon_bmi)

#chi squared test
chisq.test(x=meta_alpha[["sex"]], y=meta_alpha[["diagnosis"]])

ggplot(meta_alpha, aes(x=sex, y=diagnosis)) +
  geom_count(aes(size = after_stat(prop), group = sex)) +
  scale_size_area(max_size = 10) +
  geom_count(aes(size = after_stat(prop), group = diagnosis)) +
  scale_size_area(max_size = 10)
  scale_x_discrete(name=NULL,
                   breaks=c("female", "male"),
                   labels=c("Female", "Male")) +
  scale_y_discrete(name=NULL,
                   breaks=c("normal", "adenoma", "cancer"),
                   labels=c("Normal", "Adenoma", "Cancer")) +
  scale_size_area() +
  labs(title="There is significant variation in the likelihood that men or women will\ndevelop lesions",
       x="Body Mass Index (BMI)",
       y="Number of observed OTUs") +
  theme_classic()

