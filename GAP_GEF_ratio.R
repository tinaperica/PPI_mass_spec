library(tidyverse)

GAP_GEF_log2FC <- read_tsv("apms_log2_fold_change.txt") %>% 
  filter(gene_name %in% c("SRM1", "RNA1") & norm == "eqM") 

tag_avrg <- GAP_GEF_log2FC %>% 
  arrange(mutant, tag, gene_name) %>% 
  group_by(mutant, gene_name) %>% 
  summarise("tag_avg_log2FC" = mean(log2FC, na.rm = T), 
            "tag_avg_log2FC_sd" = sd(log2FC, na.rm = T)) %>% 
  ungroup() %>% 
  mutate("tag_avg_log2FC_sd" = ifelse(is.na(tag_avg_log2FC_sd), 0, tag_avg_log2FC_sd)) %>% 
  inner_join(., GAP_GEF_log2FC, by = c("mutant", "gene_name")) %>% 
  select(mutant, gene_name, tag_avg_log2FC, tag_avg_log2FC_sd, tag)

mut_with_both_tags <- tag_avrg %>% 
  group_by(mutant, gene_name) %>% 
  summarize("count" = n()) %>% 
  filter(count > 1) %>% pull(mutant) %>% unique()

tag_avrg <- tag_avrg %>% 
  mutate("tag"= ifelse(is.element(mutant, mut_with_both_tags), "both", tag)) %>% unique()

tag_avrg %>% ggplot(aes(x = mutant, y = tag_avg_log2FC, fill = tag)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = tag_avg_log2FC - tag_avg_log2FC_sd, ymax = tag_avg_log2FC + tag_avg_log2FC_sd), width = 0.5) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("tag averaged ln fold change") + facet_grid(~gene_name) +
  ggtitle("GAP and GEF fold change averaged for tag (sd represents a combination of N or C tag difference)")


tag_avrg <- tag_avrg %>% 
  group_by(mutant) %>% 
  summarize("gap_over_gef_FC")