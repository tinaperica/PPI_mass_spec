library(tidyverse)
library(plotly)
library(heatmaply)
(ms_input <- read_tsv("20180428_FoldChange_APMS_data.txt", col_names = T))
ms_input <- ms_input %>%
  mutate("sample" = str_c(Tag, "_", Mutant))
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
ms_data <- inner_join(ms_input, orf_gene_index, by = "ORF") %>%
  inner_join(., localization, by = "ORF")
mutants_by_fold_change_of_nuclear_proteins <- ms_data %>%
  filter(compartment == "nucleus" & Background != "NoTag") %>% 
  group_by(sample) %>%
  summarise("sum" = sum(iLog2FC)) %>%
  arrange(sum) %>% pull(sample)
ms_data %>% 
  mutate("sample" = factor(sample, mutants_by_fold_change_of_nuclear_proteins)) %>%
  filter(Background != "NoTag") %>%
  ggplot(mapping = aes(sample, iLog2FC, fill = compartment)) + 
  geom_col() + coord_flip() + ylab("Fold Change") + ggtitle("Localization of prey proteins")

### merge with interface informations
core_partners <- interface_residues %>% pull(partner) %>% unique()
interface_residues_completed <- interface_residues %>%
  select(partner, yeastresnum, mean_n_contacts) %>%
  complete(partner, nesting(yeastresnum), fill = list(mean_n_contacts = 0))
samples_ordered_by_residue_number <- ms_data %>%
  group_by(sample, Mutant) %>%
  summarise("resnum" = mean(residue)) %>%
  arrange(resnum, Mutant) %>% pull(sample)
ms_by_interface_comparison <- inner_join(ms_data, interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  filter(gene_name %in% core_partners & Background != "NoTag" & gene_name == partner) %>%
  mutate("sample" = factor(sample, samples_ordered_by_residue_number))
ggplot(ms_by_interface_comparison, aes(x = gene_name, y = sample, fill = iLog2FC, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Interfaces and interactions")

ggplot(ms_by_interface_comparison, aes(x = gene_name, y = sample, size = iPvalue, fill = iLog2FC)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Altered partner interactions per mutant") +
  scale_size("iPvalue", range = c(6, 2), breaks = c(0, 0.25, 0.5, 0.75, 1))


ms_by_interface_comparison %>% filter(iPvalue < 0.05) %>%
  ggplot(aes(x = gene_name, y = sample, size = iPvalue, fill = iLog2FC)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + ggtitle("Altered partner interactions per mutant") +
  scale_size("iPvalue", range = c(6, 2), breaks = c(0, 0.0125, 0.025, 0.0375, 0.05))
 