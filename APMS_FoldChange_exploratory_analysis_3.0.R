library(tidyverse)
library(plotly)
library(heatmaply)
library(GGally)
library(ggrepel)
library(gtools)
library(ggcorrplot)

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
#ms_data_gene_name <- inner_join(ms_input, orf_gene_index, by = "ORF") %>%
# inner_join(., localization, by = "ORF")
ms_data_all <- inner_join(ms_input, localization, by = "ORF")

### pool of real prey proteins
pool_of_all_preys <- ms_data_all %>% 
  pull(gene_name) %>% unique()
pool_of_sig_preys <- ms_data_all %>%
  filter(Background == "NoTag" & iPvalue < 0.05 & iLog2FC > 0) %>% 
  arrange(gene_name) %>% 
  pull(gene_name) %>% unique()
pool_of_sig_very_high_FC_preys <- ms_data_all %>% 
  filter(Background == "NoTag" & iPvalue < 0.05 & iLog2FC > 20)
ms_data <- ms_data_all %>% 
  filter(gene_name %in% pool_of_sig_preys & Background != "NoTag")


mutants_by_fold_change_of_nuclear_proteins <- ms_data %>%
  filter(compartment == "nucleus") %>% 
  group_by(sample) %>%
  summarise("sum" = sum(iLog2FC)) %>%
  arrange(sum) %>% pull(sample)
ms_data %>% 
  filter(iPvalue < 0.05) %>% 
  mutate("sample" = factor(sample, mutants_by_fold_change_of_nuclear_proteins)) %>%
  ggplot(mapping = aes(sample, iLog2FC, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of log2 fold change") + ggtitle("Localization of prey proteins")

avrg_loc <- ms_data %>% 
  filter(iPvalue < 0.05) %>% 
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
  geom_col() + coord_flip() + ggtitle("MeanFoldChanges", subtitle = "nuclear minus cytoplasmic prey proteins")
  
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
  complete(partner, nesting(yeastresnum), fill = list(mean_n_contacts = 1))

samples_ordered_by_residue_number <- ms_data %>%
  group_by(sample, Mutant) %>%
  summarise("resnum" = mean(residue)) %>%
  arrange(resnum, Mutant) %>% pull(sample)



ms_by_interface_comparison <- ms_data %>% 
  filter(gene_name %in% core_partners) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -10, -10, iLog2FC)) %>% 
  inner_join(., interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  select(gene_name, partner, sample, iLog2FC, mean_n_contacts,iPvalue) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample)


ms_by_interface_comparison %>%
  mutate("iLog2FC" = ifelse(iLog2FC < -10, -10, iLog2FC)) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC > 10, 10, iLog2FC)) %>% 
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, fill = iLog2FC, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Interfaces and interactions")

ms_by_interface_comparison  %>%
  mutate("sample" = factor(sample, samples_ordered_by_residue_number)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC < -10, -10, iLog2FC)) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC > 10, 10, iLog2FC)) %>% 
  ggplot(aes(x = gene_name, y = sample, size = iPvalue, fill = iLog2FC)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Altered partner interactions per mutant") +
  scale_size("iPvalue", range = c(6, 0.011), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))

#ms_by_interface_comparison %>% 
 # select(sample, iLog2FC, gene_name) %>%
#  spread(gene_name, iLog2FC) %>% 
 # select(-sample) %>% 
  #ggcorr(nbreaks = 10, label = T) + ggtitle("Pearson correlation between core partners")

#### correlations between partners
# Pearson
ms_data_spread <- ms_data %>% 
  filter(gene_name %in% core_partners) %>% 
  select(gene_name, sample, iLog2FC) %>% 
  spread(gene_name, iLog2FC)
ms_data_matrix <- as.matrix(ms_data_spread[, -1])
rownames(ms_data_matrix) <- ms_data_spread$sample
##### correlations between mutants
cormat <- cor(ms_data_matrix[,], use = "pairwise.complete.obs")
cormat %>% ggcorrplot(hc.order = TRUE,
                       outline.col = "white",
                       insig = "blank") + ggtitle("Pearson correlation")
# spearman correlation
cormat <- cor(ms_data_matrix[,], use = "pairwise.complete.obs", method = "spearman")
cormat %>% ggcorrplot(hc.order = TRUE,
                      outline.col = "white",
                      insig = "blank") + ggtitle("Kendall rank correlation")


plots <- list()
core_partner_pairs <- combn(core_partners, 2)
core_partner_pairs <- core_partner_pairs[, c(20, 22, 32, 80, 84, 99, 114)]
for (i in seq_along(core_partner_pairs[1,])) {
  partner1 <- core_partner_pairs[1, i]
  partner2 <- core_partner_pairs[2, i]
  to_plot <- ms_data %>% 
    filter(gene_name == partner1 | gene_name == partner2) %>%
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
        to_plot$meaniPvalue < 0.05 & (to_plot$partner1 > 3 | to_plot$partner2 > 3), 
        ],
      inherit.aes = F, aes(partner1, partner2, label = sample, color = Tag),
      box.padding = 0.25, point.padding = 0.25, size = 2,
      show.legend = F) + theme_bw()
  }
}
print(plots)



plots <- list()
example_mutants <- c("N_T34E", "N_T34L", "C_T34L", "N_R108L", "C_R108L", "N_H141V","C_H141V" )
#sample_pairs <- combn(unique(ms_data$sample), 2)
sample_pairs <- combn(example_mutants, 2)
for (i in seq_along(sample_pairs[1,])) {
  sample1 <- sample_pairs[1, i]
  sample2 <- sample_pairs[2, i]
  to_plot <- ms_data %>%
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
      scale_size("meaniPvalue", range = c(3, 0.01), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))
    plots[[i]] <- plots[[i]] + geom_label_repel(
      data = to_plot[to_plot$meaniPvalue < 0.01 & (abs(to_plot$sample1) > 5 | abs(to_plot$sample2) > 5), ],
      inherit.aes = F, aes(sample1, sample2, label = gene_name),
      box.padding = 0.25, point.padding = 0.25, size = 2,
      show.legend = F) + theme_bw()
  }
}
print(plots)



partners <- ms_by_interface_comparison %>% 
  select(sample, gene_name, iLog2FC) %>% 
  spread(gene_name, iLog2FC)
baits.scaled.matrix <- as.matrix(partners[, 2:9])
baits.scaled.matrix[is.na(baits.scaled.matrix)] <- 0
rownames(baits.scaled.matrix) <- partners$sample
heatmaply(baits.scaled.matrix, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85)


all_pulled_down_proteins <- ms_data %>% pull(ORF) %>% unique()



### now do the analysis separately by tag

### N-terminally tagged
### pull of real prey proteins
N_pool_of_all_preys <- ms_data_all %>% 
  filter(Tag == "N") %>% 
  pull(gene_name) %>% unique()
N_pool_of_sig_preys <- ms_data_all %>%
  filter(Tag == "N" & Background == "NoTag" & iPvalue < 0.05 & iLog2FC > 0) %>% 
  pull(gene_name) %>% unique()

ms_data_N <- ms_data_all %>% 
  filter(gene_name %in% N_pool_of_sig_preys & Tag == "N" & Background != "NoTag")

### C-terminally tagged
### pull of real prey proteins
C_pool_of_all_preys <- ms_data_all %>% 
  filter(Tag == "C") %>% 
  pull(gene_name) %>% unique()
C_pool_of_sig_preys <- ms_data_all %>%
  filter(Tag == "C" & Background == "NoTag" & iPvalue < 0.05 & iLog2FC > 0) %>% 
  pull(gene_name) %>% unique()

ms_data_C <- ms_data_all %>% 
  filter(gene_name %in% C_pool_of_sig_preys & Tag == "C" & Background != "NoTag")

##### common pool between N and C terminal tags
common_pool <- intersect(N_pool_of_sig_preys, C_pool_of_sig_preys) %>% sort()

### spread
N_ms_data_spread <- ms_data_N %>% 
  select(gene_name, sample, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -20, -20, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 20, 20, iLog2FC)) %>%
  spread(gene_name, iLog2FC)

N_ms_data_matrix <- as.matrix(N_ms_data_spread[, -1])
rownames(N_ms_data_matrix) <- N_ms_data_spread$sample
##### correlations between mutants
Ncormat_rank <- cor(t(N_ms_data_matrix[,]), use = "pairwise.complete.obs", method = "kendall")
Ncormat_rank %>% ggcorrplot(hc.order = TRUE,
           outline.col = "white",
           insig = "blank") + ggtitle("Kendall rank correlation", subtitle = "N-tagged mutants")
Ncormat <- cor(t(N_ms_data_matrix[,]), use = "pairwise.complete.obs")
Ncormat %>% ggcorrplot(hc.order = TRUE,
              outline.col = "white",
              insig = "blank") + ggtitle("Pearson correlation", "N-tagged mutants")

#### partners and mutants heatmap
#### remove columns (prey proteins) that have more than 50% NA
N_ms_data_spread_for_heatmap <- N_ms_data_spread[, -which(colMeans(is.na(N_ms_data_spread)) > 0.6)]
N_ms_data_matrix_for_heatmap <- as.matrix(N_ms_data_spread_for_heatmap[, -1])
rownames(N_ms_data_matrix_for_heatmap) <- N_ms_data_spread_for_heatmap$sample
heatmaply(N_ms_data_matrix_for_heatmap, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.2)
#### filter all interactions by iPvalue and remake heatmap
N_ms_data_spread_0.05 <- ms_data_N %>% 
  filter(iPvalue < 0.05) %>% 
  select(gene_name, sample, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -20, -20, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 20, 20, iLog2FC)) %>%
  spread(gene_name, iLog2FC)
#### remove columns that are more than 60% NA
N_ms_data_spread_for_heatmap_0.05 <- N_ms_data_spread_0.05[, -which(colMeans(is.na(N_ms_data_spread_0.05)) > 0.6)]
N_ms_data_matrix_for_heatmap_0.05 <- as.matrix(N_ms_data_spread_for_heatmap_0.05[, -1])
rownames(N_ms_data_matrix_for_heatmap_0.05) <- N_ms_data_spread_for_heatmap_0.05$sample
heatmaply(N_ms_data_matrix_for_heatmap_0.05, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.5)

C_ms_data_spread <- ms_data_C %>% 
  select(gene_name, sample, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -20, -20, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 20, 20, iLog2FC)) %>%
  spread(gene_name, iLog2FC)

C_ms_data_matrix <- as.matrix(C_ms_data_spread[, -1])
rownames(C_ms_data_matrix) <- C_ms_data_spread$sample
##### correlations between mutants
Ccormat_rank <- cor(t(C_ms_data_matrix[,]), use = "pairwise.complete.obs", method = "kendall")
Ccormat_rank %>% ggcorrplot(hc.order = TRUE,
                            outline.col = "white",
                            insig = "blank") + ggtitle("Kendall rank correlation", subtitle = "C-tagged mutants")
Ccormat <- cor(t(C_ms_data_matrix[,]), use = "pairwise.complete.obs")
Ccormat %>% ggcorrplot(hc.order = TRUE,
                 outline.col = "white",
                insig = "blank") + ggtitle("Pearson correlation", "C-tagged mutants")

#### partners and mutants heatmap
#### remove columns (prey proteins) that have more than 50% NA
C_ms_data_spread_for_heatmap <- C_ms_data_spread[, -which(colMeans(is.na(C_ms_data_spread)) > 0.8)]
C_ms_data_matrix_for_heatmap <- as.matrix(C_ms_data_spread_for_heatmap[, -1])
rownames(C_ms_data_matrix_for_heatmap) <- C_ms_data_spread_for_heatmap$sample
heatmaply(C_ms_data_matrix_for_heatmap, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.2)
#### filter all interactions by iPvalue and remake heatmap
C_ms_data_spread_0.05 <- ms_data_C %>% 
  filter(iPvalue < 0.05) %>% 
  select(gene_name, sample, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -20, -20, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 20, 20, iLog2FC)) %>%
  spread(gene_name, iLog2FC)
C_ms_data_spread_for_heatmap_0.05 <- C_ms_data_spread_0.05[, -which(colMeans(is.na(C_ms_data_spread_0.05)) > 0.6)]
C_ms_data_matrix_for_heatmap_0.05 <- as.matrix(C_ms_data_spread_for_heatmap_0.05[, -1])
rownames(C_ms_data_matrix_for_heatmap_0.05) <- C_ms_data_spread_for_heatmap_0.05$sample
heatmaply(C_ms_data_matrix_for_heatmap_0.05, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.6)


##### common pool between N and C terminal tags
common_pool <- intersect(N_pool_of_sig_preys, C_pool_of_sig_preys) %>% sort()
common_ms_data <- ms_data %>% 
  filter(gene_name %in% common_pool)


common_ms_data_spread <- common_ms_data %>% 
  select(gene_name, sample, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -20, -20, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 20, 20, iLog2FC)) %>%
  spread(gene_name, iLog2FC)

common_ms_data_matrix <- as.matrix(common_ms_data_spread[, -1])
rownames(common_ms_data_matrix) <- common_ms_data_spread$sample
##### correlations between mutants
common_cormat_rank <- cor(t(common_ms_data_matrix[,]), use = "pairwise.complete.obs", method = "kendall")
common_cormat_rank %>% ggcorrplot(hc.order = TRUE,
                            outline.col = "white",
                            insig = "blank") + ggtitle("Kendall rank correlation", subtitle = "common preys")
common_cormat <- cor(t(common_ms_data_matrix[,]), use = "pairwise.complete.obs")
common_cormat %>% ggcorrplot(hc.order = TRUE,
                       outline.col = "white",
                       insig = "blank") + ggtitle("Pearson correlation", "common preys")

#### partners and mutants heatmap
#### remove columns (prey proteins) that have more than 50% NA
common_ms_data_spread_for_heatmap <- common_ms_data_spread[, -which(colMeans(is.na(common_ms_data_spread)) > 0.8)]
common_ms_data_matrix_for_heatmap <- as.matrix(common_ms_data_spread_for_heatmap[, -1])
rownames(common_ms_data_matrix_for_heatmap) <- common_ms_data_spread_for_heatmap$sample
heatmaply(common_ms_data_matrix_for_heatmap, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.2)



#### wodak complexes
ms_data_N_wodak <- inner_join(ms_data_N, wodak, by = "ORF") %>% 
  select(sample, gene_name, complex, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -30, -30, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 30, 30, iLog2FC)) %>%
  group_by(sample, complex) %>% 
  summarise("mean_per_complex_FC" = mean(iLog2FC, na.rm = T)) %>% 
  ungroup() %>% 
  spread(complex, mean_per_complex_FC)
ms_data_N_wodak_for_heatmap <- ms_data_N_wodak[, -which(colMeans(is.na(ms_data_N_wodak)) > 0.9)]
ms_data_N_wodak_matrix <- as.matrix(ms_data_N_wodak_for_heatmap[, -1])
rownames(ms_data_N_wodak_matrix) <- ms_data_N_wodak$sample
ms_data_C_wodak <- inner_join(ms_data_C, wodak, by = "ORF") %>% 
  select(sample, gene_name, complex, iLog2FC) %>% 
  mutate("iLog2FC" = ifelse(iLog2FC < -30, -30, iLog2FC)) %>%
  mutate("iLog2FC" = ifelse(iLog2FC > 30, 30, iLog2FC)) %>%
  group_by(sample, complex) %>% 
  summarise("mean_per_complex_FC" = mean(iLog2FC, na.rm = T)) %>% 
  ungroup() %>% 
  spread(complex, mean_per_complex_FC)
ms_data_C_wodak_for_heatmap <- ms_data_C_wodak[, -which(colMeans(is.na(ms_data_C_wodak)) > 0.9)]
ms_data_C_wodak_matrix <- as.matrix(ms_data_C_wodak_for_heatmap[, -1])
rownames(ms_data_C_wodak_matrix) <- ms_data_C_wodak$sample


heatmaply(ms_data_N_wodak_matrix, margins = c(160, 70), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.8)
heatmaply(ms_data_C_wodak_matrix, margins = c(160, 70), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, cexCol = 0.8)
