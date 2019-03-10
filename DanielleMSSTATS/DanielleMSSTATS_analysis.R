library(tidyverse)
library(ggcorrplot)
library(RColorBrewer)
library(pcaMethods)
library(ggrepel)
library(ggfortify)
# z_score <- function(value, mean, sd) {
#   z_score <- (value - mean) / sd
#   return(z_score)
# }
### input data
N_PPIs <- read_tsv("DanielleMSSTATS/SAINT_all_N-term_preys_5percBFDR.txt",col_names = T)
C_PPIs <- read_tsv("DanielleMSSTATS/SAINT_all_C-term_preys_5percBFDR.txt", col_names = T)
### use the msstats results based only on SAINT identified preys (that means use only preys confirmed with SAINT)
### log2 fold change files
N_eqM_log2 <- read_tsv("DanielleMSSTATS/N-ter_PPI_eqM_msstats_results.txt", col_names = T) %>% ### equalized median
  mutate("tag" = "N", "norm" = "eqM")
C_eqM_log2 <- read_tsv("DanielleMSSTATS/C-ter_PPI_eqM_msstats_results.txt", col_names = T) %>% 
  mutate("tag" = "C", "norm" = "eqM")
N_gs_log2 <- read_tsv("DanielleMSSTATS/N-ter_PPI_gs_msstats_results.txt", col_names = T) %>% ### global standard
  mutate("tag" = "N", "norm" = "gs")
C_gs_log2 <- read_tsv("DanielleMSSTATS/C-ter_PPI_gs_msstats_results.txt", col_names = T) %>% 
  mutate("tag" = "C", "norm" = "gs")
### parse and prepare the data
ms_log2 <- bind_rows(N_eqM_log2, C_eqM_log2, N_gs_log2, C_gs_log2)
ms_log2 <- ms_log2 %>% 
  mutate('residue' = as.numeric(substring(Label, 2, (nchar(Label)-1)))) %>% 
  select("PreyORF" = Protein, "mutant" = Label, log2FC, adj.pvalue, residue, tag, norm) %>% 
  mutate("sample" = str_c(tag, mutant, sep = "_")) %>% 
  mutate("log2FC" = -1 * log2FC)    ##### convert fold change from wt/mut to mut/wt

### abundance files
N_eqM_abundance <- read_tsv("DanielleMSSTATS/Gsp1-N-eqM_results-mss-groupQuant.txt", col_names = T) %>% 
  gather(key = mutant, value = abundance, -Protein) %>% 
  mutate("tag" = "N", "norm" = "eqM")
C_eqM_abundance <- read_tsv("DanielleMSSTATS/Gsp1-C-eqM_results-mss-groupQuant.txt", col_names = T) %>% 
  gather(key = mutant, value = abundance, -Protein) %>% 
  mutate("tag" = "C", "norm" = "eqM")
N_gs_abundance <- read_tsv("DanielleMSSTATS/Gsp1-N-gs_results-mss-groupQuant.txt", col_names = T) %>% 
  gather(key = mutant, value = abundance, -Protein) %>% 
  mutate("tag" = "N", "norm" = "gs")
C_gs_abundance <- read_tsv("DanielleMSSTATS/Gsp1-C-gs_results-mss-groupQuant.txt", col_names = T) %>% 
  gather(key = mutant, value = abundance, -Protein) %>% 
  mutate("tag" = "C", "norm" = "gs")
ms_abundance <- bind_rows(N_eqM_abundance, N_gs_abundance, C_eqM_abundance, C_gs_abundance) %>% 
  rename("PreyORF" = Protein) %>% 
  mutate('residue' = as.numeric(substring(mutant, 2, (nchar(mutant)-1)))) %>% 
  mutate("sample" = str_c(tag, mutant, sep = "_")) %>% 
  group_by(norm, tag) %>% 
  mutate("z_score" = scale(abundance, center = T, scale = T)) %>%   ### this scaling is the same as z-score
  ungroup()
### input auxiliary files
#### index to match ORFs to gene names
sgd_orf <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F)
orf_gene_index <- sgd_orf %>% select("ORF" = X1, "gene_name" = X2) %>% unique()
### wodak complexes
wodak <- read_tsv("wodak_complex.txt", col_names = T, 
                  col_type = list(gene_name = col_skip()))
### GO-slims file
GO_slims <- read_tsv("20180216_go_slim_mapping.tab.txt", col_names = F, 
                     col_types=list(X1 = 'c', X2 = 'c', X3 = "_", X4 = 'c', X5 = 'c', X6 = '_', X7 = '_'))
names(GO_slims) <- c("ORF", "gene_name", "GO_type", "GO")

### interfaces
interfaces <- read_tsv("DanielleMSSTATS/SASA_interfaces.txt", col_names = T)
core_partners <- interfaces %>% pull(partner) %>% unique()
interface_residues_completed <- interfaces %>%
  filter(interface == "core") %>% 
  select(partner, yeastresnum, deltarASA) %>%
  complete(partner, nesting(yeastresnum), fill = list(deltarASA = 0)) %>% 
  arrange(yeastresnum) %>% 
  group_by(partner) %>% 
  mutate(deltaASA_z_score = scale(deltarASA, center = T, scale = T)) %>% 
  ungroup()

all_interface_residues_completed <- interfaces %>%
  group_by(partner, yeastresnum) %>% 
  mutate("deltarASA" = mean(deltarASA)) %>% 
  ungroup() %>% 
  select(partner, yeastresnum, deltarASA) %>%
  complete(partner, nesting(yeastresnum), fill = list(deltarASA = 0)) %>% 
  arrange(yeastresnum) %>% 
  group_by(partner) %>% 
  mutate(deltaASA_z_score = scale(deltarASA, center = T, scale = T)) %>% 
  ungroup()

ms_log2 <- ms_log2 %>% 
  inner_join(., orf_gene_index, by = c("PreyORF" = "ORF"))
write_tsv(ms_log2, "apms_log2_fold_change.txt")

ms_log2 %>% 
  ggplot(aes(x = log2FC)) + 
  geom_histogram() + 
  facet_grid(tag~norm) +
  ggtitle("Distribution of log2FC values for the two normalization method")

ms_abundance <- ms_abundance %>% 
  inner_join(., orf_gene_index, by = c("PreyORF" = "ORF"))

write_tsv(ms_abundance, "APMS_Gsp1_abundance.txt")
ms_abundance %>% 
  ggplot(aes(x = abundance)) +
  geom_histogram() + 
  facet_grid(tag ~ norm) + 
  ggtitle("Distribution of prey abundances for the two normalization methods")
ggsave("DanielleMSSTATS/Abundance_distributions.pdf", width = 7, height = 7)
samples_ordered_by_residue_number <- ms_log2 %>%
  group_by(sample, mutant) %>%
  summarise("resnum" = mean(residue)) %>%
  arrange(resnum, mutant) %>% pull(sample)
ms_by_interface_comparison <- ms_log2 %>% 
  filter(gene_name %in% core_partners) %>% 
  inner_join(., all_interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  select(norm, tag, gene_name, partner, sample, log2FC, deltarASA, adj.pvalue) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample) %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number))

ms_by_interface_comparison %>%
  ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = deltarASA)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  facet_grid(~norm) +
  ggtitle("Interfaces and interactions (mutant/WT log2)")
ggsave("DanielleMSSTATS/Interfaces_and_interactions.pdf", width = 10, height = 7)
core_partners <- append(core_partners, c("GSP1", "MOG1"))
ms_log2 %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>% 
  arrange(sample) %>% 
  filter(gene_name %in% core_partners) %>% 
  ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = adj.pvalue)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  scale_size("adj.pvalue", range = c(6, 0.1), breaks = c(0, 0.001, 0.0375, 0.05, 0.1, 0.25, 0.5)) +
  facet_wrap(~norm) +
  ggtitle("Interactions with core partners (mut/wt log2)")
##### main observation from this is that it's better to use eqM because at least the same mutant with different tags look similar
ggsave("DanielleMSSTATS/Interactions_with_core_partners_log2FC.pdf", height = 8, width = 11)


ms_by_interface_comparison %>%
  filter(norm == "eqM") %>% 
  ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = deltarASA)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  labs(fill = 'log2\nfold change') +
  xlab("prey protein") + ylab ("Gsp1 mutant") +
  ggtitle("Interfaces and interactions")
ggsave("DanielleMSSTATS/Interfaces_and_interactions_eqM.pdf", width = 7, height = 8)
ms_by_interface_comparison %>% 
  filter(norm == "eqM") %>%
  mutate(adj.pvalue = ifelse(adj.pvalue > 0.05, 0.11, adj.pvalue)) %>% 
  ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = adj.pvalue)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  scale_size("p-value", range = c(6, 0.1), breaks = c(0, 0.001, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15)) +
  labs(fill = 'log2\nfold change') +
  xlab("prey protein") + ylab ("Gsp1 mutant") +
  ggtitle("Interactions with core partners (mut/wt log2)")
##### main observation from this is that it's better to use eqM because at least the same mutant with different tags look similar
ggsave("DanielleMSSTATS/Interactions_with_core_partners_log2FC_eqM.pdf", height = 8, width = 7)


# int <- ms_by_interface_comparison %>%
#   filter(norm == "eqM") %>% 
#   mutate("plot" = "interface") %>%
#   select(sample, partner, log2FC, deltarASA, plot) %>% 
#   rename("value" = deltarASA)
# msstats <- ms_by_interface_comparison %>% 
#   filter(norm == "eqM") %>% 
#   mutate("plot" = "msstats") %>% 
#   select(sample, partner, log2FC, adj.pvalue, plot) %>%
#   rename("value" = adj.pvalue)
#   bind_rows(., int)
# 
# ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = deltarASA)) +
#   geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
#   facet_grid(~plot) +
#   ggplot(aes(x = gene_name, y = sample, fill = log2FC, size = adj.pvalue)) +
#   geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
#   scale_size("adj.pvalue", range = c(6, 0.1), breaks = c(0, 0.001, 0.0375, 0.05, 0.1, 0.25, 0.5)) +

##### main observation from this is that it's better to use eqM because at least the same mutant with different tags look similar
ggsave("DanielleMSSTATS/Interactions_with_core_partners_log2FC.pdf", height = 8, width = 11)



#### look at variation in abundance with overlap in interfaces
### simple null hypothesis is that if there is overlap in the interface abundances will be
## inverseley correlated
### first calculate percent abundance for each partner (which is abundance / max abundance (for that prey for all the mutants))
ms_abundance <- ms_abundance %>%
  arrange(PreyORF) %>% 
  group_by(PreyORF, tag, norm) %>% 
  mutate("max_abundance" = max(abundance, na.rm = T)) %>% 
  ungroup() %>% 
  mutate("percent_abundance" = abundance/max_abundance)

# abundance_and_interfaces <- ms_abundance %>% 
#   filter(gene_name %in% core_partners) %>% 
#   select(tag, norm, sample, gene_name, percent_abundance) %>% 
#   inner_join(., interface_residues_completed, by = c("gene_name" = "partner")) %>% 
#   select(tag, norm, sample, gene_name, percent_abundance, deltarASA, yeastresnum)
prey_core_proteins <- ms_abundance %>% 
  filter(gene_name %in% core_partners) %>% 
  #select(tag, norm, sample, gene_name, percent_abundance) %>% 
  inner_join(., interface_residues_completed, by = c("gene_name" = "partner")) %>% 
  pull(gene_name) %>% unique()
#prey_core_protein_pairs <- combn(x = prey_core_proteins, m = 2)
prey_core_protein_pairs <- expand.grid(prey_core_proteins, prey_core_proteins)
pairwise_abundance <- tibble("tag" = 'c', "norm" = 'c', "sample" = 'c', "gene_name.x" = 'c',
                          "percent_abundance.x" = double(), "z_score.x" = double(), "gene_name.y" = 'c', 
                          "percent_abundance.y" = double(), "z_score.y" = double())
pairwise_interface <- tibble("gene_name.x" = 'c', "yeastresnum" = integer(), "deltarASA.x" = double(),
                             "deltaASA_z_score.x" = double(), 
                             "gene_name.y" = 'c', "deltarASA.y" = double(),  "deltaASA_z_score.y" = double()
                             )
for (i in seq_along(prey_core_protein_pairs[,1])) {
  prot1 <- prey_core_protein_pairs[i, 1]
  prot2 <- prey_core_protein_pairs[i, 2]
  temp1 <- ms_abundance %>% 
    filter(gene_name == prot1) %>% 
    select(tag, norm, sample, residue, gene_name, percent_abundance, z_score)  
  temp2 <- ms_abundance %>% 
    filter(gene_name == prot2) %>% 
    select(tag, norm, sample, residue, gene_name, percent_abundance, z_score)
  merged <- inner_join(temp1, temp2, by = c("tag", "norm", "sample", "residue"))
  pairwise_abundance <- bind_rows(pairwise_abundance, merged)
  temp1 <- all_interface_residues_completed %>% 
    filter(partner == prot1) %>% rename("gene_name" = partner)
  temp2 <- all_interface_residues_completed %>% 
    filter(partner == prot2) %>% rename("gene_name" = partner)
  merged <- inner_join(temp1, temp2, by = "yeastresnum")
  pairwise_interface <- bind_rows(pairwise_interface, merged)
}

plot_abundance_interfaces <- function(nor) {
  to_plot_abundance <- pairwise_abundance %>% 
    filter(norm == nor)
  to_plot_abundance %>% ggplot(aes(x = percent_abundance.x, y = percent_abundance.y, color = tag)) +
    geom_abline(intercept = 0, slope = 1, size = 0.1) +
    geom_point(alpha = 0.5, size = 1) + facet_grid(gene_name.y ~ gene_name.x) +
    xlab("percent") + ylab("percent") +
    theme_bw()  +
    geom_point(data = pairwise_interface, aes(x = deltarASA.x, y = deltarASA.y), alpha = 0.5, size = 1, color = "#e79f00")
  ggsave(str_c("DanielleMSSTATS/Abundance_and_deltarASA_for_partner_pairs_", nor, ".pdf"), height = 10, width = 14)
}
plot_abundance_interfaces("gs")
plot_abundance_interfaces("eqM")


abundance_and_SASA_for_mut_resi <- pairwise_abundance %>% 
  select(tag, norm, sample, gene_name.x, percent_abundance.x, percent_abundance.y, gene_name.y, residue) %>% 
  inner_join(., pairwise_interface[, c(1:3, 5:6)], by = c("residue" = "yeastresnum", "gene_name.x", "gene_name.y"))

pairwise_abundance_gap_gef_yrb1 <- pairwise_abundance %>% 
  filter(gene_name.x %in% c("RNA1", "SRM1", "YRB1") & gene_name.y %in% c("RNA1", "SRM1", "YRB1"))

pairwise_interface_gap_gef_yrb1 <- pairwise_interface %>% 
  filter(gene_name.x %in% c("RNA1", "SRM1", "YRB1") & gene_name.y %in% c("RNA1", "SRM1", "YRB1"))
  
abundance_and_SASA_for_mut_resi_gap_gef_yrb1 <- pairwise_abundance_gap_gef_yrb1 %>% 
  select(tag, norm, sample, gene_name.x, percent_abundance.x, percent_abundance.y, gene_name.y, residue) %>% 
  inner_join(., pairwise_interface[, c(1:3, 5:6)], by = c("residue" = "yeastresnum", "gene_name.x", "gene_name.y"))
plot_abundance_interfaces_with_mut <- function(nor) {
  to_plot_abundance <- pairwise_abundance_gap_gef_yrb1 %>% 
    filter(norm == nor & ! is.na(residue))
  wt_abundance_to_plot <- pairwise_abundance_gap_gef_yrb1 %>% 
    filter(norm == nor & is.na(residue))
  to_plot_abundance_SASA_mut <- abundance_and_SASA_for_mut_resi_gap_gef_yrb1 %>% 
    filter(norm == nor) %>% 
    filter(deltarASA.x > 0.3 | deltarASA.y > 0.3)
  to_plot_abundance %>% ggplot(aes(x = percent_abundance.x, y = percent_abundance.y, color = as.character(residue))) +
    geom_abline(intercept = 0, slope = 1, size = 0.1) +
    labs(x = "percent", y = "percent", color = "residue") +
    geom_point(alpha = 0.4, size = 1) + facet_grid(gene_name.y ~ gene_name.x) +
    xlab("percent") + ylab("percent") +
    theme_bw()  +
    geom_point(data = pairwise_interface_gap_gef_yrb1, aes(x = deltarASA.x, y = deltarASA.y), alpha = 0.3, size = 1, color = "#e79f00") +
    geom_point(data = to_plot_abundance_SASA_mut, aes(x = percent_abundance.x, y = percent_abundance.y), alpha = 1, size = 1.4) +
    geom_point(data = to_plot_abundance_SASA_mut, aes(x = deltarASA.x, y = deltarASA.y, color = as.character(residue)), alpha = 1, size = 1.4) +
    geom_point(data = wt_abundance_to_plot, aes(x = percent_abundance.x, y = percent_abundance.y), color = "black", size = 1.4)
  
  ggsave(str_c("DanielleMSSTATS/Abundance_and_deltarASA_for_partner_pairs_with_mut_gap_gef_yrb1_", nor, ".pdf"), height = 6, width = 9)
}
plot_abundance_interfaces_with_mut("gs")
plot_abundance_interfaces_with_mut("eqM")


samples_ordered_by_residue_number_plus_wt <- ms_abundance %>%
  group_by(sample, mutant) %>%
  summarise("resnum" = mean(residue)) %>%
  arrange(resnum, mutant) %>% pull(sample)
ms_abundance %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number_plus_wt)) %>% 
  arrange(sample) %>% 
  filter(gene_name %in% core_partners) %>% 
  ggplot(aes(x = gene_name, y = sample, fill = z_score, size = z_score)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  facet_wrap(~norm) +
  ggtitle("Abundance of core partners (as z-score))")
##### main observation from this is that it's better to use eqM because at least the same mutant with different tags look similar
ggsave("DanielleMSSTATS/Core_partners_abundance.pdf", height = 8, width = 9)


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
get_reordered_cormat <- function(data) {
  cormat <- round(cor(data[, -1], use = "pairwise.complete.obs", method = "spearman"), 2)
  cormat <- reorder_cormat(cormat)
  return(cormat)
}
get_order_by_corr <- function(data) {
  cormat <- get_reordered_cormat(data)
  ordered <- rownames(cormat)
}
make_log2FC_heatmap <- function(data, nor, t, width_inches, outfile_prefix) {
  if (t == "both") {
    ms_data_tag <- data %>%
      filter(norm == nor) %>% 
      mutate("mutant" = sample) %>% 
      select(mutant, log2FC, gene_name)
  } else {
    ms_data_tag <- data %>%
      filter(tag == t & norm == nor) %>% 
      select(mutant, log2FC, gene_name) %>% 
      arrange(gene_name, mutant)
  }
  spread_preys <- ms_data_tag %>% 
    spread(gene_name, log2FC)
  ordered_preys <- get_order_by_corr(spread_preys)
  spread_mutants <- ms_data_tag %>% 
    spread(mutant, log2FC)
  ordered_mutants <- get_order_by_corr(spread_mutants)
  ms_data_tag %>% 
    mutate("mutant" = factor(mutant, ordered_mutants)) %>%
    mutate("gene_name" = factor(gene_name, ordered_preys)) %>% 
    arrange(mutant, gene_name) %>% 
    ggplot(aes(x = gene_name, y = mutant, fill = log2FC)) +
    geom_tile() + scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(str_c(outfile_prefix, "_", t, "_", nor, ".pdf"), width = width_inches)
  
}

make_log2FC_heatmap(ms_log2, "eqM", "N", 40, "DanielleMSSTATS/All_interactions_log2FC")
make_log2FC_heatmap(ms_log2, "gs", "N", 40, "DanielleMSSTATS/All_interactions_log2FC")
make_log2FC_heatmap(ms_log2, "eqM", "C", 25, "DanielleMSSTATS/All_interactions_log2FC")
make_log2FC_heatmap(ms_log2, "gs", "C", 25, "DanielleMSSTATS/All_interactions_log2FC")
make_log2FC_heatmap(ms_log2, "gs", "both", 45, "DanielleMSSTATS/All_interactions_log2FC")


ms_data_complex <- ms_log2 %>% 
  inner_join(., wodak, by = c("PreyORF" = "ORF")) %>% 
  arrange(cluster) %>% 
  group_by(sample, cluster, mutant, tag, norm) %>% 
  summarize("log2FC" = mean(log2FC, na.rm = T)) %>% 
  rename("gene_name" = cluster) %>% 
  ungroup()
make_log2FC_heatmap(ms_data_complex, "eqM", "C", 10, "DanielleMSSTATS/Interactions_by_complexes_log2FC")
make_log2FC_heatmap(ms_data_complex, "eqM", "N", 10, "DanielleMSSTATS/Interactions_by_complexes_log2FC")
make_log2FC_heatmap(ms_data_complex, "gs", "C", 10, "DanielleMSSTATS/Interactions_by_complexes_log2FC")
make_log2FC_heatmap(ms_data_complex, "gs", "N", 10, "DanielleMSSTATS/Interactions_by_complexes_log2FC")

#### cluster and correlate mutants and prey genes based on abundance
mutant_corr_by_abundance <- function(nor, t) {
  prey_abundance_mat <- ms_abundance %>% 
    filter(tag == t & norm == nor) %>% 
    select(mutant, gene_name, z_score) %>% 
    spread(gene_name, z_score)
  #prey_cormat <- reorder_cormat(round(cor(prey_abundance_mat[, -1], use = "pairwise.complete.obs"), 2))
  prey_cormat <- get_reordered_cormat(prey_abundance_mat)
  prey_cormat %>% ggcorrplot(
                        outline.col = "white", colors = c("white", "#e79f00", "#0072B2"),
                        insig = "blank") + ggtitle("Pearson correlation based on prey abundance")
  ggsave(str_c("DanielleMSSTATS/Prey_abundace_zscore_corr_", t, "_", nor, ".pdf"), width = 40, height = 40)
  mutant_abundance_mat <- ms_abundance %>% 
    filter(tag == t & norm == nor) %>% 
    select(mutant, gene_name, z_score) %>% 
    spread(mutant, z_score)
  #mutant_cormat <- reorder_cormat(round(cor(mutant_abundance_mat[, -1], use = "pairwise.complete.obs"), 2))
  mutant_cormat <- get_reordered_cormat(mutant_abundance_mat)
  #mutant_cormat <- as_tibble(mutant_cormat_mat, rownames = "mutant")
  #mutant_order <- rownames(mutant_cormat_mat) 
  mutant_cormat %>% 
    ggcorrplot(
    outline.col = "white", colors = c("white", "#e79f00", "#0072B2"),
    insig = "blank") + ggtitle(str_c("Spearman correlation based on prey abundance - ", t, " tag" ))
  # mutant_cormat %>% 
  #   mutate("mutant" = factor(mutant, mutant_order)) %>% 
  #   gather(key = mutant1, value = corr, -mutant) %>% 
  #   ggplot(aes(x = mutant, y = mutant1, fill = corr)) +
  #   geom_tile() + scale_fill_gradientn(colors = brewer.pal(9, "OrRd")) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   ylab("mutant")
  ggsave(str_c("DanielleMSSTATS/Mutants_abundace_zscore_corr_", t, "_", nor, ".pdf"), width = 7, height = 7)
}
mutant_corr_by_abundance("eqM", "N")
mutant_corr_by_abundance("eqM", "C")
mutant_corr_by_abundance("gs", "N")
mutant_corr_by_abundance("gs", "C")

get_dd <- function(data) {
  cormat <- cor(data[,-1], use = "pairwise.complete.obs")
  dd <- as.dist((1 - cormat)/2, upper = T, diag = T)
  dd[is.na(dd)] <- 0
  return(dd)
}

abundance_heatmap <- function(nor, t) {
  temp <- ms_abundance %>% 
    filter(tag == t & norm == nor) %>% 
    select(mutant, gene_name, z_score)
  mutants <- temp %>% 
    pull(mutant) %>% unique()
  preys <- temp %>% 
    pull(gene_name) %>% unique()
  mut_spread <- temp %>% 
    spread(mutant, z_score)
  dd <- get_dd(mut_spread)
  ordered_mutants <- mutants[hclust(dd)$order]
  prey_spread <- temp %>% 
    spread(gene_name, z_score)
  dd <- get_dd(prey_spread)
  ordered_preys <- preys[hclust(dd)$order]
  temp %>% 
    mutate("mutant" = factor(mutant, ordered_mutants),
           "gene_name" = factor(gene_name, ordered_preys)) %>% 
    ggplot(aes(x= gene_name, y = mutant, fill = z_score)) +
    geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "PiYG")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(str_c("DanielleMSSTATS/Mutants_vs_prey_abundace_zscore_heatmap_", t, "_", nor, ".pdf"), width = 32)
  
}
abundance_heatmap("eqM", "N")
abundance_heatmap("eqM", "C")
abundance_heatmap("gs", "N")
abundance_heatmap("gs", "C")

### z-scores between the two normalization methods
eqM <- ms_abundance %>% 
  filter(norm == "eqM") %>% 
  select(-norm)
gs <- ms_abundance %>% 
  filter(norm == "gs") %>% 
  select(-norm)
merge_eqM_gs <- inner_join(eqM, gs, by = names(eqM)[c(1:2, 4:6, 8)])
merge_eqM_gs %>% 
  ggplot(aes(x = abundance.x, abundance.y)) + geom_point(color = "gray") +
  geom_text(aes(label = ifelse((abundance.x/abundance.y < 0.91 | abundance.y/abundance.x < 0.91), str_c(mutant, " ", gene_name), '')), hjust = -0.1, vjust = 1) +
  xlab("abundance in eqM") + ylab("abundance in gs") + geom_abline(slope = 1, intercept = 0)

merge_eqM_gs %>% 
  ggplot(aes(x = z_score.x, z_score.y)) + geom_point(color = "gray") +
  xlab("z-score in eqM") + ylab("z-score in gs") + geom_abline(slope = 1, intercept = 0)

N_abundance <- ms_abundance %>% 
  filter(tag == "N") %>% 
  select(mutant, norm, z_score, gene_name) %>% 
  mutate("z_score" = ifelse(is.na(z_score), 0, z_score))
C_abundance <- ms_abundance %>% 
  filter(tag == "C") %>% 
  select(mutant, norm, z_score, gene_name) %>% 
  mutate("z_score" = ifelse(is.na(z_score), 0, z_score))
both_tag_mutants <- N_abundance %>% 
  inner_join(., C_abundance, by = c("mutant", "norm", "gene_name")) %>% 
  filter(! is.na(z_score.x) & ! is.na(z_score.y) )

both_tag_mutants <- both_tag_mutants %>% 
  group_by(norm) %>% 
  mutate("predicted" = predict(lm(z_score.y ~ z_score.x))) %>% 
  mutate("residual" = residuals(lm(z_score.y ~ z_score.x))) %>% 
  mutate("abs_residual" = abs(residual)) %>% 
  arrange(desc(abs_residual))
top_preys_diff_between_tags <- both_tag_mutants %>% 
  filter(abs_residual > 2) %>% 
  pull(gene_name) %>% unique()
to_color <- both_tag_mutants %>% 
  filter(gene_name %in% union(c("SRM1", "YRB1"), top_preys_diff_between_tags))
# both_tag_mutants %>% 
#   ggplot(aes(x = z_score.x, y = z_score.y)) + 
#   geom_point() + 
#   geom_point(data = to_color, aes(color = gene_name)) +
#   #geom_text(aes(label = ifelse(abs_residual > 2, gene_name, ''), hjust = -0.1, vjust = 1)) +
#   geom_line(aes(x = z_score.x, y = predicted)) +
#   xlab("z-score N tag") + ylab("z score C tag") + 
#   xlim(c(-4, 4)) + ylim(c(-4, 4)) +
#   ggtitle("Abundance for N vs C tag (mutants with both tags)") +
#   facet_grid(~norm)

#select_mut <- c("WT", "G80A", "K143W", "K132H")
both_tag_mutants %>% 
  #filter(mutant %in% select_mut) %>% 
  ggplot(aes(x = z_score.x, y = z_score.y)) +  #, shape = mutant)) + 
  geom_point() + 
  geom_point(data = to_color[to_color$mutant %in% select_mut,], aes(color = gene_name)) +
  #geom_text(aes(label = ifelse(abs_residual > 2, gene_name, ''), hjust = -0.1, vjust = 1)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("z-score N tag") + ylab("z score C tag") + 
  xlim(c(-4, 4)) + ylim(c(-4, 4)) +
  ggtitle("Abundance for N vs C tag - mutants with both tags only") +
  facet_grid(~norm)


top_preys_diff_between_tags <- both_tag_mutants %>% 
  filter(abs_residual > 1.5) %>% 
  pull(gene_name) %>% unique()
to_color <- both_tag_mutants %>% 
  filter(gene_name %in% union(c("SRM1", "YRB1"),  top_preys_diff_between_tags))
both_tag_mutants %>% 
  filter(mutant == "WT") %>% 
  ggplot(aes(x = z_score.x, y = z_score.y)) + 
  geom_point() + 
  geom_point(data = to_color[to_color$mutant == "WT",], aes(color = gene_name)) +
  #geom_text(aes(label = ifelse(abs_residual > 2, gene_name, ''), hjust = -0.1, vjust = 1)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("z-score N tag") + ylab("z score C tag") + 
  xlim(c(-4, 4)) + ylim(c(-4, 4)) +
  ggtitle("Abundance for N vs C tag WT only") +
  facet_grid(~norm)
### principal component analysis
ms_abundance_for_pca <- ms_abundance %>% 
  filter(norm == "eqM") %>%
  select(mutant, tag, z_score, gene_name) %>% 
  spread(gene_name, z_score)
ms_abundance_for_pca <- data.frame(ms_abundance_for_pca)
ms_abundance_for_pca[is.na(ms_abundance_for_pca)] <- 0

autoplot(prcomp(ms_abundance_for_pca[, -c(1, 2)]), data = ms_abundance_for_pca, colour = 'tag')
ggsave("DanielleMSSTATS/PCA_all_preys_both_tags.pdf")

ms_abundance_for_pca <- ms_abundance %>% 
  filter(! gene_name %in% top_preys_diff_between_tags & norm == "eqM") %>% 
  select(mutant, tag, z_score, gene_name) %>% 
  spread(gene_name, z_score)
ms_abundance_for_pca <- data.frame(ms_abundance_for_pca)
ms_abundance_for_pca[is.na(ms_abundance_for_pca)] <- 0
ggsave("DanielleMSSTATS/PCA_all_preys_both_tags_top_diff_preys_removed.pdf")



ms_abundance_for_pca <- ms_abundance %>% 
  filter(norm == "eqM") %>%
  select(mutant, tag, z_score, gene_name) %>% 
  spread(gene_name, z_score)
#ms_abundance_for_pca <- data.frame(ms_abundance_for_pca)
ms_abundance_for_pca[is.na(ms_abundance_for_pca)] <- 0

autoplot(prcomp(ms_abundance_for_pca[, -c(1, 2)]), data = ms_abundance_for_pca, colour = 'tag')
ggsave("DanielleMSSTATS/PCA_all_preys_both_tags.pdf")





ms_abundance_for_pca <- ms_abundance %>% 
  filter(norm == "eqM" & tag == "N") %>%
  #filter(gene_name %in% core_partners) %>% 
  filter(gene_name != "GSP1") %>% 
  select(mutant, z_score, gene_name) %>% 
  spread(gene_name, z_score)
ms_abundance_for_pca_both_tags <- ms_abundance %>% 
  filter(norm == "eqM") %>%
  #filter(gene_name %in% core_partners) %>% 
  filter(gene_name != "GSP1") %>% 
  select("mutant" = sample, z_score, gene_name) %>% 
  spread(gene_name, z_score)
pca_wrap <- function(data, pcn = 2, tag) {
  #### first remove columns with more than half NA
  #data <- data[, -which(colMeans(is.na(data)) > 0.5)]  # remove columns that are more than half NA
  data_preped <- pcaMethods::prep(data[, -1], center = FALSE, scale = "none")
  pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  scores <- as_tibble(scale(scores(pca))) %>% 
    add_column(., "mutant" = data$mutant)
  loadings <- as_tibble(scale(loadings(pca)))
  loadings <- add_column(loadings, "variables" = colnames(data[-1]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings, "cum_var" = R2cum(pca)))
}

principal_components <- str_c("PC", 1:7)
pca <- pca_wrap(ms_abundance_for_pca_both_tags, pcn = length(principal_components))
#core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
pca$loadings <- pca$loadings %>% 
  mutate("PreyGene" = str_sub(variables, -4))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = pca$scores, aes_string("PC1", pc)) + 
    #geom_label(aes(label = mutant), label.size = 0.1) +
    geom_point() +
    geom_label_repel(aes(label = mutant)) +
    theme_classic() +
    #geom_point(data = pca$loadings, aes_string("PC1", pc, color = "PreyGene"), 
              #size = 3, alpha = 0.5)
    geom_point(data = pca$loadings, aes_string("PC1", pc), 
      size = 1, alpha = 0.1)
}
pdf("DanielleMSSTATS/Both_tags_PCA.pdf")
print(plots)
dev.off()




ms_abundance_for_pca_both_tags <- ms_abundance %>% 
  filter(norm == "eqM") %>%
  #filter(gene_name %in% core_partners) %>% 
  filter(gene_name != "GSP1") %>% 
  select(mutant, tag, z_score, gene_name) %>% 
  spread(gene_name, z_score)
pca_wrap <- function(data, pcn = 2) {
  #### first remove columns with more than half NA
  #data <- data[, -which(colMeans(is.na(data)) > 0.5)]  # remove columns that are more than half NA
  data_preped <- pcaMethods::prep(data[, c(-1, -2)], center = FALSE, scale = "none")
  #pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  data_preped[is.na(data_preped)] <- 0
  pca <- prcomp(data_preped, scale. = T)
  plot(pca)
  scores <- as_tibble(scale(pca$x)) %>% 
    add_column(., "mutant" = data$mutant) %>% 
    add_column(., "tag" = data$tag)
  loadings <- as_tibble(scale(pca$rotation))
  loadings <- add_column(loadings, "variables" = colnames(data[, c(-1, -2)]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings))
}

principal_components <- str_c("PC", 1:10)
pca <- pca_wrap(ms_abundance_for_pca_both_tags, pcn = length(principal_components))
#core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
#pca$loadings <- pca$loadings %>% 
 # mutate("PreyGene" = str_sub(variables, -4))
df <- pca$loadings[, 1:39]
df[df < 0.5 & df > -0.5] <- NA
df <- df %>% 
  mutate("variables" = pca$loadings$variables)
#df <- cbind(df, "variables" = as.character(pca$loadings$variables))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = pca$scores, aes_string("PC2", pc, color = "tag")) + 
    #geom_label(aes(label = mutant), label.size = 0.1) +
    geom_point() +
    geom_label_repel(aes(label = mutant)) +
    theme_classic() +
    geom_point(data = pca$loadings, aes_string("PC1", pc), 
    size = 3, alpha = 0.1)
    #geom_point(data = pca$loadings, aes_string("PC1", pc), 
     #          size = 1, alpha = 0.1)
}
pdf("DanielleMSSTATS/Both_tags_PCA_PC2.pdf")
print(plots)
dev.off()



ms_abundance_for_pca_both_tags <- ms_abundance %>% 
  filter(norm == "eqM") %>%
  filter(gene_name %in% core_partners) %>% 
  filter(gene_name != "GSP1") %>% 
  select(mutant, tag, z_score, gene_name) %>% 
  spread(gene_name, z_score)
pca_wrap <- function(data, pcn = 2) {
  #### first remove columns with more than half NA
  #data <- data[, -which(colMeans(is.na(data)) > 0.5)]  # remove columns that are more than half NA
  data_preped <- pcaMethods::prep(data[, c(-1, -2)], center = FALSE, scale = "none")
  #pca <- pcaMethods::pca(data_preped, nPcs = pcn)
  data_preped[is.na(data_preped)] <- 0
  pca <- prcomp(data_preped, scale. = T)
  plot(pca)
  scores <- as_tibble(scale(pca$x)) %>% 
    add_column(., "mutant" = data$mutant) %>% 
    add_column(., "tag" = data$tag)
  loadings <- as_tibble(scale(pca$rotation))
  loadings <- add_column(loadings, "variables" = colnames(data[, c(-1, -2)]))
  return(list("all" = pca, "scores" = scores, "loadings" = loadings))
}

principal_components <- str_c("PC", 1:10)
pca <- pca_wrap(ms_abundance_for_pca_both_tags, pcn = length(principal_components))
#core_pca <- pca_wrap(core_regulator_data, pcn = length(principal_components))
#pca$loadings <- pca$loadings %>% 
# mutate("PreyGene" = str_sub(variables, -4))
plots <- list()
for (pc in principal_components[2:length(principal_components)]) {
  plots[[pc]] <-  ggplot(data = pca$scores, aes_string("PC1", pc, color = "tag")) + 
    #geom_label(aes(label = mutant), label.size = 0.1) +
    geom_point() +
    geom_label_repel(aes(label = mutant)) +
    theme_classic() +
    geom_point(data = pca$loadings, aes_string("PC1", pc, color = "variables"), 
               size = 3, alpha = 0.5)
  #geom_point(data = pca$loadings, aes_string("PC1", pc), 
  #          size = 1, alpha = 0.1)
}
pdf("DanielleMSSTATS/Both_tags_PCA_core_partners.pdf")
print(plots)
dev.off()


