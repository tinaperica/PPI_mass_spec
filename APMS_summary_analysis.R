library(tidyverse)
apms <- read_delim("20180322_clean_APMS_summary.txt", delim = "\t", col_names = T)
#samples_from_Jan <- apms %>% filter(`Month of collection` == "January") %>%
 # pull(`Sample ID`)
#apms %>% filter(`Sample ID` %in% samples_from_Jan) %>%
 # ggplot(aes(x = `Total peptide IDs`, color = `Month of collection`)) + geom_density()
ggsave("total_peptides_per_month.pdf")
apms %>% ggplot(aes(x = `Total peptide IDs`, color = as.character(pulldown))) + geom_density()
ggsave("total_peptides_per_pulldown_exp.pdf")
apms %>% filter(`Total peptide IDs` > 1000) %>%
  ggplot(aes(y = `Condition`, x = `Total peptide IDs`)) + geom_point()
ggsave("total_peptide_IDs_per_sample.pdf", height = 20)

apms %>% filter(`Total peptide IDs` > 1000) %>%
  ggplot(aes(y = `Condition`, x = `Total peptide IDs`, color = as.character(pulldown))) + geom_point()
ggsave("total_peptide_IDs_per_sample_col_by_pulldown.pdf", height = 10)


apms <- apms %>% filter(`Gsp1 intensity` > 1e11) %>%
  mutate("mean_Gsp1_intensity" = mean(`Gsp1 intensity`)) %>%
  arrange(mean_Gsp1_intensity)
order_by_mean_Gsp1_intensity <- apms %>% pull(Condition) %>% unique()
apms %>% 
  mutate(sample = factor(Condition, levels = order_by_mean_Gsp1_intensity)) %>%
  ggplot(aes(y = `Condition`, x = `Gsp1 intensity`, color = as.character(pulldown))) + geom_point()
ggsave("Gsp1_intensity_per_sample_col_by_pulldown_March_only_cleaned_data.pdf", height = 10)


apms %>% filter(`Total peptide IDs` > 1000) %>%
  ggplot(aes(y = `Condition`, x = `Total peptide IDs`, color = `Tag terminus`)) + geom_point()


apms <- apms %>% filter(`Total peptide IDs` > 2000) %>%
  group_by(Condition) %>%
  mutate("peptide IDs range" = max(`Total peptide IDs`) - min(`Total peptide IDs`),
         "peptide IDs SD" = sd(`Total peptide IDs`)) %>%
  arrange(`peptide IDs SD`)
order_by_peptide_SD <- apms %>% pull(`Condition`) %>% unique()
apms %>%
  mutate(sample = factor(Condition, levels = order_by_peptide_SD)) %>%
  ggplot(aes(y = sample, x = `peptide IDs SD`)) + geom_point()
ggsave("per_sample_peptide_IDs_stdev_March_only_cleaned_data.pdf")

apms <- apms %>% filter(`Gsp1 intensity` > 0) %>%
  group_by(Condition) %>%
  mutate("Gsp1 intensity SD" = sd(`Gsp1 intensity`)) %>%
  arrange(`Gsp1 intensity SD`)
order_by_Gsp1_intensity <- apms %>% pull(`Condition`) %>% unique()
apms %>%
  mutate(sample = factor(Condition, levels = order_by_Gsp1_intensity)) %>%
  ggplot(aes(y = sample, x = `Gsp1 intensity SD`)) + geom_point()
ggsave("per_sample_Gsp1_intensity_stdev_March_only_cleaned_data.pdf")

apms %>% filter(`Gsp1 intensity` > 0 & `Total peptide IDs` > 0) %>%
  ggplot(aes(x = `Total peptide IDs`, y = `Gsp1 intensity`)) + geom_point()

apms %>% ggplot(aes(x = `Total peptide IDs`, y = Run)) + geom_point()

                