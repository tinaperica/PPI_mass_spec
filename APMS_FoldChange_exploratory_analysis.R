library(tidyverse)
library(plotly)
library(heatmaply)
library(GGally)
library(ggrepel)
library(gtools)

(ms_input <- read_tsv("20180502_FoldChange_APMS_data.txt", col_names = T))
ms_input <- ms_input %>%
  mutate("sample" = str_c(Tag, "_", Mutant))
unique(ms_input$sample)
#### index to match ORFs to gene names
sgd_orf <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F)
orf_gene_index <- sgd_orf %>% select("ORF" = X1, "gene_name" = X2) %>% unique()

### wodak complexes
wodak <- read_delim("wodak_complex.txt", delim = "\t", col_names = T, 
                    col_type = list(gene_name = col_skip()))
wodak <- wodak %>%
  rename("complex" = cluster)

##### get interface information for the mutants
interface_residues <- read_tsv("contacts_and_interfaces.txt")

#### extract the residue
ms_input <- ms_input %>%
  mutate("residue" = as.numeric(substring(Mutant, 2, (nchar(Mutant)-1))))
### get localization information on prey proteins
### is it in the nucleus, cytoplasm, or both
GO_slims <- read_delim("20180216_go_slim_mapping.tab.txt", delim = "\t", col_names = F)
localization <- GO_slims %>% 
  filter(X4 == "C" & (X5 == "nucleus" | X5 == "cytoplasm")) %>%
  select(c(1,5)) %>%
  rename("ORF" = X1, "compartment" = X5)
nuc_cyto_genes <- localization %>% 
  group_by(ORF) %>%
  summarize("count" = n()) %>%
  filter(count == 2) %>% pull(ORF)
uniq_loc_genes <- localization %>%
  filter(! ORF %in% nuc_cyto_genes)
ambig_loc_genes <- tibble("ORF" = nuc_cyto_genes, compartment = "nuc/cyto")
localization <- bind_rows(uniq_loc_genes, ambig_loc_genes)
currated_loc <- read_tsv("localization_manual_currated.txt")
currated_loc <- inner_join(currated_loc, orf_gene_index, by = c("PreyGeneName" = "gene_name")) %>%
  select(c(3, 2))
localization <- localization %>%
  filter(! ORF %in% currated_loc$ORF) %>%
  bind_rows(., currated_loc)


##### add gene names and localization
ms_data_gene_name <- inner_join(ms_input, orf_gene_index, by = "ORF") %>%
  inner_join(., localization, by = "ORF")
### filter out by pvalue (but keep all the imputed hits)
ms_data <- ms_data_gene_name %>%
  filter(Background != "NoTag") %>%
  filter(iPvalue < 0.05)
mutants_by_fold_change_of_nuclear_proteins <- ms_data %>%
  filter(compartment == "nucleus") %>% 
  group_by(sample) %>%
  summarise("sum" = sum(iLog2FC)) %>%
  arrange(sum) %>% pull(sample)
ms_data %>% 
  mutate("sample" = factor(sample, mutants_by_fold_change_of_nuclear_proteins)) %>%
  ggplot(mapping = aes(sample, iLog2FC, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of log2 fold change") + ggtitle("Localization of prey proteins")

avrg_loc <- ms_data %>% 
  filter(compartment == "nucleus" | compartment == "cytoplasm") %>%
  group_by(sample, compartment) %>%
  summarize("average log2 fold change" = mean(iLog2FC)) %>%
  ungroup() %>%
  spread(compartment, `average log2 fold change`) %>%
  mutate("nuc minus cyto log2" = nucleus - cytoplasm)
samples_by_avrg_nuc_loc <- avrg_loc %>%
  arrange(nucleus) %>% 
  pull(sample)
avrg_loc %>%
  mutate("sample" = factor(sample, samples_by_avrg_nuc_loc)) %>%
  gather(key = "compartment", value = "log2", - `nuc minus cyto log2`, -sample) %>%
  ggplot(aes(sample, log2, fill = compartment)) +
  geom_col() + coord_flip()
samples_by_nuc_ratio <- avrg_loc %>%
  arrange(`nuc minus cyto log2`) %>%
  pull(sample)
avrg_loc %>% 
  mutate("sample" = factor(sample, samples_by_nuc_ratio)) %>%
  ggplot(aes(sample, `nuc minus cyto log2`)) + geom_col() + coord_flip()

### merge with interface informations
core_partners <- interface_residues %>% pull(partner) %>% unique()
interface_residues_completed <- interface_residues %>%
  select(partner, yeastresnum, mean_n_contacts) %>%
  complete(partner, nesting(yeastresnum), fill = list(mean_n_contacts = 0.5))

samples_ordered_by_residue_number <- ms_data %>%
  group_by(sample, Mutant) %>%
  summarise("resnum" = mean(residue)) %>%
  arrange(resnum, Mutant) %>% pull(sample)



ms_by_interface_comparison <- ms_data %>% 
  filter(gene_name %in% core_partners & Background != "NoTag") %>% 
  #complete(gene_name, nesting(sample, residue), fill = list(iLog2FC = -10, iPvalue = NA)) %>% 
  inner_join(., interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  select(gene_name, partner, sample, iLog2FC, mean_n_contacts,iPvalue) %>% 
  #complete(gene_name, partner, nesting(sample), fill = list(iLog2FC = -10, mean_n_contacts = 0.5)) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample)


ms_by_interface_comparison %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, fill = iLog2FC, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Interfaces and interactions")

ms_by_interface_comparison  %>%
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, size = iPvalue, fill = iLog2FC)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Altered partner interactions per mutant") +
  scale_size("iPvalue", range = c(6, 0.011), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))

ms_by_interface_comparison %>% 
  select(sample, iLog2FC, gene_name) %>%
  spread(gene_name, iLog2FC) %>% 
  select(-sample) %>% 
  ggcorr(nbreaks = 10, label = T) + ggtitle("Pearson correlation between core partners")
plots <- list()
core_partner_pairs <- combn(core_partners, 2)
for (i in seq_along(core_partner_pairs[1,])) {
  partner1 <- core_partner_pairs[1, i]
  partner2 <- core_partner_pairs[2, i]
  to_plot <- ms_data_gene_name %>% 
      filter(Background != "NoTag" & (gene_name == partner1 | gene_name == partner2)) %>%
      group_by(sample) %>%
      mutate("meaniPvalue" = mean(iPvalue)) %>%
      ungroup() %>%
      select(sample, gene_name, iLog2FC, meaniPvalue, Tag) %>%
      spread(gene_name, iLog2FC)
  if (ncol(to_plot) == 5) {
    names(to_plot)[4:5] <- c("partner1", "partner2")
    plots[[i]] <- ggplot(to_plot, aes(partner1, partner2, size = meaniPvalue, color = Tag)) + 
      geom_point() + ggtitle(str_c(partner1, partner2, sep = " ")) +
      xlab(partner1) + ylab(partner2) + xlim(-10, 10) + ylim(-10, 10) +
      scale_size("meaniPvalue", range = c(6, 1), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05, 0.1))
    plots[[i]] <- plots[[i]] + geom_label_repel(
      data = to_plot[ 
        (to_plot$partner1 > 5 | to_plot$partner2 > 5) |
           ( ( to_plot$partner1 > 2 | to_plot$partner2 > 2 ) & 
              to_plot$partner1/to_plot$partner2 < 3 & to_plot$partner1/to_plot$partner2 > 0.6 ) , 
          ],
      inherit.aes = F, aes(partner1, partner2, label = sample, color = Tag),
      box.padding = 0.25, point.padding = 0.25, size = 2,
      show.legend = F) + theme_bw()
  }
}
print(plots)



plots <- list()
example_mutants <- c("N_T34E", "N_T34A", "N_T34L", "C_T34L", "C_T34Q", "N_R108L", "C_R108L", "N_H141V","C_H141V" )
#sample_pairs <- combn(unique(ms_data$sample), 2)
sample_pairs <- combn(example_mutants, 2)
for (i in seq_along(sample_pairs[1,])) {
  sample1 <- sample_pairs[1, i]
  sample2 <- sample_pairs[2, i]
  to_plot <- ms_data_gene_name %>%
    filter(Background != "NoTag" & (sample == sample1 | sample == sample2)) %>%
    group_by(gene_name) %>%
    mutate("meaniPvalue" = mean(iPvalue)) %>%
    ungroup() %>%
    select(sample, gene_name, iLog2FC, meaniPvalue, compartment) %>%
    spread(sample, iLog2FC)
  if(ncol(to_plot) == 5) {
    names(to_plot)[4:5] <- c("sample1", "sample2")
    plots[[i]] <- ggplot(to_plot, aes(sample1, sample2, size = meaniPvalue, color = compartment)) +
      geom_point() + ggtitle(str_c(sample1, sample2, sep = " ")) +
      xlab(sample1) + ylab(sample2) + xlim(-13, 13) + ylim(-13, 13) +
      scale_size("meaniPvalue", range = c(3, 0.01), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))
    plots[[i]] <- plots[[i]] + geom_label_repel(
      data = to_plot[to_plot$sample1 > 7 | to_plot$sample2 > 7, ],
      inherit.aes = F, aes(sample1, sample2, label = gene_name),
      box.padding = 0.25, point.padding = 0.25, size = 2,
      show.legend = F) + theme_bw()
  }
}
print(plots)


ms_data %>% 
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 4, layout.exp = 3) +
  ggtitle("Correlation between mutants based on all interactions", 
  subtitle = "(FDR < 0.05, normalized by tagged WT)")

ms_data %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 4, layout.exp = 3) +
  ggtitle("Correlation between mutants based on all interactions", 
          subtitle = "(FDR < 0.05, normalized by tagged WT)")


ms_by_interface_comparison %>% 
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 3, layout.exp = 3) +
  ggtitle("Correlation between mutant based on only core partner interactions", 
          subtitle = "(FDR < 0.05, normalized by tagged WT)")

ms_by_interface_comparison %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 3, layout.exp = 3) +
  ggtitle("Correlation between mutant based on only core partner interactions", 
          subtitle = "(FDR < 0.05, normalized by tagged WT)")


partners <- ms_by_interface_comparison %>% 
  select(sample, gene_name, iLog2FC) %>% 
  spread(gene_name, iLog2FC)
baits.scaled.matrix <- as.matrix(partners[, 2:13])
baits.scaled.matrix[is.na(baits.scaled.matrix)] <- 0
rownames(baits.scaled.matrix) <- partners$sample
heatmaply(baits.scaled.matrix, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85)


#### remove all the prey genes that have negative iLog2FC against the NoTag
genes_enriched_in_background <- ms_data_gene_name %>%
  filter(Background == "NoTag" & iLog2FC < 0) %>%
  pull(gene_name) %>% unique()
ms_data_no_background <- ms_data_gene_name %>%
  filter(iPvalue < 0.05) %>%
  filter(! gene_name %in% genes_enriched_in_background) %>%
  filter(Background != "NoTag")

ms_data_no_background %>%
  #mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 4, layout.exp = 3) +
  ggtitle("Correlation between mutants based on all interactions", 
          subtitle = "(FDR < 0.05, normalized by tagged WT, all background genes removed)")


plots <- list()
example_mutants <- c("N_T34E", "N_T34A", "N_T34L", "C_T34L", "C_T34Q", "N_R108L", "C_R108L", "N_H141V","C_H141V" )
sample_pairs <- combn(example_mutants, 2)
#sample_pairs <- combn(unique(ms_data_no_background$sample), 2)
for (i in seq_along(sample_pairs[1,])) {
  sample1 <- sample_pairs[1, i]
  sample2 <- sample_pairs[2, i]
  to_plot <- ms_data_no_background %>% 
    filter(sample == sample1 | sample == sample2) %>%
    group_by(gene_name) %>%
    mutate("meaniPvalue" = mean(iPvalue)) %>%
    ungroup() %>%
    select(sample, gene_name, iLog2FC, meaniPvalue, compartment) %>%
    spread(sample, iLog2FC)
  if(ncol(to_plot) == 5) {
    names(to_plot)[4:5] <- c("sample1", "sample2")
    plots[[i]] <- ggplot(to_plot, aes(sample1, sample2, size = meaniPvalue, color = compartment)) + 
      geom_point() + ggtitle(str_c(sample1, sample2, sep = " ")) +
      xlab(sample1) + ylab(sample2) + xlim(-13, 13) + ylim(-13, 13) +
      scale_size("meaniPvalue", range = c(6, 1), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))
    plots[[i]] <- plots[[i]] + geom_label_repel(
      data = to_plot,
      inherit.aes = F, aes(sample1, sample2, label = gene_name),
      box.padding = 0.25, point.padding = 0.25, size = 2,
      show.legend = F) + theme_bw()
  }
}
print(plots)


######
# for each mutant prey pair check that the prey is significant for the mut/background

ms_data_wt_control <- ms_data_gene_name %>%
  filter(Background != "NoTag")
ms_data_background_control <- ms_data_gene_name %>%
  filter(Background == "NoTag") %>%
  select(sample, gene_name, "BackiLog2FC" = iLog2FC, "BackiPvalue" = iPvalue)

ms_data_wt_background <- inner_join(ms_data_wt_control, ms_data_background_control, by = c("sample", "gene_name"))
ms_data_more_than_neg_control <- ms_data_wt_background %>%
  filter(BackiPvalue < 0.05 & BackiLog2FC > 0)
ms_data_less_than_or_same_neg_control <- ms_data_wt_background %>%
  filter(BackiLog2FC <  0 | BackiPvalue > 0.05)

ms_data_more_than_neg_control %>%
  select(sample, iLog2FC, gene_name) %>%
  spread(sample, iLog2FC) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 4, layout.exp = 3) +
  ggtitle("Correlation between mutants based on all interactions", 
          subtitle = "(normalized by tagged WT, but all background genes removed)")



#### test "scenarios"
x <- c("+", "-", "not sig")
all_scenarios <- as_tibble(
  permutations(n = 3, r = 3, v = x, repeats.allowed = T)
  ) %>%
  rename("mut/wt" = V1, "mut/back" = V2, "wt/back" = V3)
(feasable_depleted_scenarios <- all_scenarios %>%
  slice(c(2, 3, 5, 8)))
(feasable_enriched_scenarios <- all_scenarios %>% 
    slice(c(10, 13:16, 18))
    )
(feasable_same_as_wt_scenarios <- all_scenarios %>% 
    slice(c(19, 23))
)
