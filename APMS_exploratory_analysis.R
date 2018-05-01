library(tidyverse)
library(plotly)
library(heatmaply)
(ms_data <- read_delim("20180402_apms_results.txt", delim = "\t", col_names = T))

#### index to match ORFs to gene names
sgd_orf <- read_delim("orf_gene_GO_sgd_annotation.txt", delim = "\t", col_names = F)
orf_gene_index <- sgd_orf %>% select("ORF" = X1, "gene_name" = X2) %>% unique()

### wodak complexes
wodak <- read_delim("wodak_complex.txt", delim = "\t", col_names = T, 
                    col_type = list(gene_name = col_skip()))
wodak <- wodak %>%
  rename("complex" = cluster)

##### get interface information for the mutants
interface_res_num <- read_tsv("interface_residues.txt")
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

tags <- paste(c("c", "n"), collapse = "|")
tags <- paste0("(", tags, ")")
ms_data <- ms_data %>% 
  mutate("Mut" = gsub(Bait, pattern = "Gsp1_", replacement = "")) %>%
  mutate(Mut = sub(tags, "\\1|", Mut)) %>%
  separate(Mut, into = c("tag", "mutant"))

##### add gene names and localization
ms_data <- inner_join(ms_data, orf_gene_index, by = c("Prey" = "ORF")) %>%
  select(Bait, tag, mutant, "PreyORF" = Prey, "PreyGeneName" = gene_name, Intensity, IntensitySum, AvgIntensity, NumReplicates, ctrlIntensity, FoldChange, BFDR)

ms_data <- inner_join(ms_data, localization, by = c("PreyORF" = "ORF")) %>%
  select(c(1:5, 13, 11:12, 6:10))

n_wt_ms_data <- ms_data %>%
  filter(Bait == "Gsp1_nWT") %>%
  select(PreyORF, "WTFoldChange" = FoldChange, "WTAvgIntensity" = AvgIntensity)
n_ms_data <- ms_data %>%
  filter(tag == "n")
n_ms_data <- inner_join(n_ms_data, n_wt_ms_data, by = "PreyORF")
c_wt_ms_data <- ms_data %>%
  filter(Bait == "Gsp1_cWT") %>%
  select(PreyORF, "WTFoldChange" = FoldChange, "WTAvgIntensity" = AvgIntensity)
c_ms_data <- ms_data %>%
  filter(tag == "c")
c_ms_data <- inner_join(c_ms_data, c_wt_ms_data, by = "PreyORF")

ms_data <- bind_rows(n_ms_data, c_ms_data) %>%
  mutate("RelFoldChange" = FoldChange / WTFoldChange, "RelAvgIntensity" = AvgIntensity / WTAvgIntensity)
### note RelAvgIntensity and RelFoldChange are basically the same number

#### add a residue number and an interface column
ms_data <- ms_data %>%
  mutate("residue" = as.numeric(substring(mutant, 2, (nchar(mutant)-1)))) %>% 
  arrange(PreyGeneName, FoldChange)
ms_data <- inner_join(ms_data, interface_res_num, by = c("residue" = "res_num"))


write_delim(ms_data, "Gsp1_PPIs_20180402_gene_names.txt", delim = "\t")

##### analyse complexes (use Wodak complex annotations)
ms_data_complex <- inner_join(ms_data, wodak, by = c("PreyORF" = "ORF")) %>%
  select(1, 4:6, 19:20, 8, 16, 2:3, 7, 9:15, 17:18) %>%
  arrange(complex, PreyORF, Bait)

####### heatmaps
### small heatmap with just the core partners
partners_to_check <- c("YRB1", "MOG1", "SRM1", "RNA1", "PSE1", "MSN5", "KAP95", 
                       "LOS1", "KAP104", "CRM1", "CSE1", "YRB2", "SRP1", "NTF2", "NUP1", "NUP60")
### baits as row numbers
partners <- ms_data %>% 
  filter(PreyGeneName %in% partners_to_check) %>%
  #filter(BFDR < 0.05) %>%
  select(Bait, PreyGeneName, AvgIntensity) %>% 
  spread(PreyGeneName, as.numeric(AvgIntensity))
baits.scaled.matrix <- as.matrix(partners[, 2:14])
rownames(baits.scaled.matrix) <- partners$Bait
heatmaply(baits.scaled.matrix, margins = c(60, 110), scale = "none",
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85)

#### interfaces as row numbers
partners <- ms_data %>%
  filter(PreyGeneName %in% partners_to_check) %>%
  select(Bait, mutant, interface, PreyGeneName, RelAvgIntensity)
matrix <- partners %>% spread(PreyGeneName, as.numeric(RelAvgIntensity))
scaled.matrix <- scale(matrix[,3:15])
rownames(scaled.matrix) <- matrix[,2]
heatmaply(matrix, margins = c(60, 110),
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, krow = 10)



partners <- ms_data %>%
  filter(PreyGeneName %in% partners_to_check & 
           (PreyGeneName == "YRB1" | PreyGeneName == "YRB2")) %>%
  select(Bait, mutant, interface, PreyGeneName, RelAvgIntensity)
matrix <- partners %>% spread(PreyGeneName, as.numeric(RelAvgIntensity))
scaled.matrix <- scale(matrix[,4])
rownames(scaled.matrix) <- matrix[,1]
heatmaply(matrix, margins = c(60, 110),
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85, krow = 10)



all.mat <- ms_data %>%
  filter(BFDR < 0.01) %>%
  select(Bait, PreyGeneName, RelAvgIntensity) %>%
  complete(Bait, PreyGeneName) %>%
  replace_na(list(RelAvgIntensity = 0)) %>%
  select(Bait, PreyGeneName, RelAvgIntensity) %>%
  spread(PreyGeneName, as.numeric(RelAvgIntensity))

baits.scaled.matrix <- scale(all.mat[-1])
rownames(baits.scaled.matrix) <- all.mat$Bait
heatmaply(baits.scaled.matrix, margins = c(60, 110),
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.85)

baits_by_prey_count <- ms_data %>%
  filter(BFDR < 0.01) %>%
  select(Bait, PreyORF) %>%
  group_by(Bait) %>%
  summarise("count" = n()) %>%
  arrange(count) %>% pull(Bait) 
ms_data %>%
  filter(BFDR < 0.01) %>%
  mutate("Bait" = factor(Bait, baits_by_prey_count)) %>%
  select(Bait, PreyGeneName, "protein localization" = compartment) %>% 
  ggplot(by.compartment, mapping = aes(Bait, fill = `protein localization`)) + 
  geom_bar() + coord_flip() + ylab("Number of proteins") + ggtitle("Localization of prey proteins - 0.01 BFDR")
ggsave("localization_of_pray_proteins_0.01.pdf")

baits_by_avg_intensity <- ms_data %>%
  filter(BFDR < 0.05) %>%
  select(Bait, RelAvgIntensity) %>%
  group_by(Bait) %>%
  summarise("IntensitySum" = sum(RelAvgIntensity)) %>%
  arrange(IntensitySum) %>%
  pull(Bait) %>% unique()
ms_data %>%
  filter(BFDR < 0.05) %>%
  mutate("Bait" = factor(Bait, baits_by_avg_intensity)) %>%
  select(Bait, compartment, RelAvgIntensity) %>%
  group_by(Bait, compartment) %>%
  ggplot(aes(Bait, RelAvgIntensity, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of AvgIntensity")
ggsave("localization_of_prey_proteins_RelAvgIntensitySum_0.05.pdf")


baits_by_RelFoldChange <- ms_data %>%
  #filter(BFDR < 0.05) %>%
  select(Bait, RelFoldChange) %>%
  group_by(Bait) %>%
  summarise("SumRelFoldChange" = sum(RelFoldChange)) %>%
  arrange(MeanRelFoldChange) %>%
  pull(Bait) %>% unique()
ms_data %>%
  #filter(BFDR < 0.05) %>%
  mutate("Bait" = factor(Bait, baits_by_RelFoldChange)) %>%
  select(Bait, compartment, RelFoldChange) %>%
  group_by(Bait, compartment) %>%
  ggplot(aes(Bait, RelFoldChange, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of FoldChange rel to WT flag")
ggsave("localization_of_prey_proteins_RelFoldChange.pdf")





ms_data_complex_0.5FDR <- ms_data_complex %>% filter(BFDR < 0.05)
per.comp.mat <- ms_data_complex %>%
  #filter(BFDR < 0.05) %>%
  select(Bait, complex, RelFoldChange) %>%
  complete(Bait, complex) %>%
  replace_na(list(RelFoldChange = 0)) %>%
  group_by(Bait, complex) %>%
  summarise("meanRelFoldChange" = mean(RelFoldChange), 
            "sumRelFoldChange" = sum(RelFoldChange)) %>%
  select(Bait, complex, sumRelFoldChange) %>%
  spread(complex, as.numeric(sumRelFoldChange))

per.comp.baits.scaled.matrix <- scale(per.comp.mat[-1])
rownames(per.comp.baits.scaled.matrix) <- per.comp.mat$Bait
heatmaply(per.comp.baits.scaled.matrix, margins = c(215, 105),
          scale_fill_gradient_fun = scale_fill_gradient2(), cexRow = 0.7, cexCol = 0.8)
             




baits_by_fold_change <- ms_data %>%
  filter(BFDR < 0.05) %>%
  select(Bait, RelFoldChange) %>%
  group_by(Bait) %>%
  summarise("RelFoldChangeSum" = sum(RelFoldChange)) %>%
  arrange(RelFoldChangeSum) %>%
  pull(Bait) %>% unique()
ms_data %>%
  filter(BFDR < 0.05) %>%
  mutate("Bait" = factor(Bait, baits_by_fold_change)) %>%
  select(Bait, compartment, RelFoldChange) %>%
  group_by(Bait, compartment) %>%
  ggplot(aes(Bait, RelFoldChange, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of Relative Fold Change")
ggsave("currated_localization_of_prey_proteins_RelFoldChangeSum_0.05.pdf")

baits_by_fold_change <- ms_data %>%
  filter(BFDR < 0.05) %>%
  select(Bait, RelFoldChange) %>%
  group_by(Bait) %>%
  summarise("RelFoldChangeSum" = sum(RelFoldChange)) %>%
  arrange(RelFoldChangeSum) %>%
  pull(Bait) %>% unique()
ms_data %>%
  filter(BFDR < 0.05) %>%
  mutate("Bait" = factor(Bait, baits_by_fold_change)) %>%
  select(Bait, compartment, RelFoldChange) %>%
  group_by(Bait, compartment) %>%
  ggplot(aes(Bait, RelFoldChange, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of Relative Fold Change")
ggsave("currated_localization_of_prey_proteins_RelFoldChangeSum_0.05.pdf")

baits_by_avg_intensity <- ms_data %>%
  filter(BFDR < 0.05) %>%
  select(Bait, AvgIntensity) %>%
  group_by(Bait) %>%
  summarise("AvgIntensitySum" = sum(AvgIntensity)) %>%
  arrange(AvgIntensitySum) %>%
  pull(Bait) %>% unique()

ms_data %>%
  filter(BFDR < 0.05) %>%
  mutate("Bait" = factor(Bait, baits_by_avg_intensity)) %>%
  select(Bait, compartment, AvgIntensity) %>%
  group_by(Bait, compartment) %>%
  ggplot(aes(Bait, AvgIntensity, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Sum of Average Intensity")
ggsave("currated_localization_of_prey_proteins_AvgIntensity_0.05.pdf")




