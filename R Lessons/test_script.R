#module 1
library(tidyverse)

r_version <- R.version$version.string

read_csv("https://raw.githubusercontent.com/riffomonas/data/master/comma-survey/comma-survey.csv") %>%
  rename(data=`How would you write the following sentence?`) %>%
  mutate(data=recode(data,
                     `Some experts say it's important to drink milk, but the data are inconclusive.` = "Plural",
                     `Some experts say it's important to drink milk, but the data is inconclusive.` = "Singular")
  ) %>%
  count(data) %>%
  drop_na() %>%
  mutate(percentage = 100 * n/sum(n)) %>%
  ggplot(aes(x=data, y=percentage, fill=data)) +
  geom_col(show.legend=FALSE) +
  labs(x=NULL,
       y="Percentage of respondents",
       title="Is the word 'data' plural or singular?",
       subtitle=r_version) +
  theme_classic()
