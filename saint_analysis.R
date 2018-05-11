library(tidyverse)
library(plotly)
library(heatmaply)
library(GGally)
library(ggrepel)

#### index to match ORFs to gene names
sgd_orf <- read_tsv("orf_gene_GO_sgd_annotation.txt", col_names = F)
orf_gene_index <- sgd_orf %>% select("ORF" = X1, "gene_name" = X2) %>% unique()

### wodak complexes
wodak <- read_tsv("wodak_complex.txt", col_names = T, 
                    col_type = list(gene_name = col_skip()))
wodak <- wodak %>%
  rename("complex" = cluster)

##### get interface information for the mutants
interface_residues <- read_tsv("contacts_and_interfaces.txt")

### GO and localization
GO_slims <- read_delim("20180216_go_slim_mapping.tab.txt", delim = "\t", col_names = F)
localization <- GO_slims %>% 
  filter(X4 == "C" & (X5 == "nucleus" | X5 == "cytoplasm")) %>%
  select(c(1,5)) %>%
  rename("ORF" = X1, "compartment" = X5)
####
saint <- read_tsv("saint_intensity.txt", col_names = T)

saint <- saint %>%
  mutate("residue" = as.numeric(substring(Mutant, 2, (nchar(Mutant)-1)))) %>% 
  mutate("sample" = str_c(Tag, "_", Mutant)) %>% 
  inner_join(orf_gene_index, by = c("PreyORF" = "ORF"))

saint <- saint %>% 
  mutate("scaled" = scale(AvgP)) 

saint %>% ggplot(aes(x = scaled)) + geom_density()
saint %>% ggplot(aes(x = AvgP)) + geom_density()

samples_by_residue_number <- saint %>% 
  group_by(sample, Mutant) %>% 
  summarise("resnum" = mean(residue)) %>% 
  arrange(resnum, Mutant) %>% 
  pull(sample)

core_partners <- interface_residues %>% pull(partner) %>% unique()
interface_residues_completed <- interface_residues %>%
  select(partner, yeastresnum, mean_n_contacts) %>%
  complete(partner, nesting(yeastresnum), fill = list(mean_n_contacts = 1))

saint_and_interfaces <- inner_join(saint, interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  select(gene_name, partner, sample, scaled, AvgP, mean_n_contacts) %>% 
  complete(gene_name, nesting(sample), fill = list(scaled = -1.5)) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample)
  
saint_and_interfaces %>% 
  mutate("sample" = factor(sample, samples_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, fill = scaled, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  ggtitle("Interfaces and interactions - SAINT intensity - scaled probability")

saint_and_interfaces %>% 
  select(gene_name, partner, sample, AvgP, mean_n_contacts, scaled) %>% 
  complete(gene_name, nesting(sample), fill = list(AvgP = -1.5)) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample)

saint_and_interfaces %>% 
  mutate("sample" = factor(sample, samples_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, fill = AvgP, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  ggtitle("Interfaces and interactions - SAINT intensity - probability")


saint_real <- saint %>% 
  mutate("AvgP" = ifelse(BFDR > 0.1, 0, AvgP)) %>%
  mutate("scaled" = scale(AvgP))

saint_real %>% ggplot(aes(x = AvgP, fill = Tag)) + geom_histogram() +
  facet_grid(~Tag)

saint_and_interfaces_sig <- inner_join(saint_real, interface_residues_completed, by = c("residue" = "yeastresnum")) %>%
  select(gene_name, partner, sample, scaled, AvgP, mean_n_contacts, BFDR) %>%
  complete(gene_name, nesting(sample), fill = list(scaled = -1.5)) %>% 
  filter(gene_name %in% core_partners & gene_name == partner) %>% 
  arrange(partner, sample)

saint_and_interfaces_sig %>% 
  mutate("sample" = factor(sample, samples_by_residue_number)) %>%
  ggplot(aes(x = gene_name, y = sample, fill = scaled, size = mean_n_contacts)) +
  geom_point(shape = 21, stroke = 0.1) + scale_fill_gradient2() + 
  ggtitle("Interfaces and interactions - SAINT intensity - scaled probability (BFDR < 0.1)")


saint_and_interfaces %>% select(sample, gene_name, AvgP) %>% 
  spread(sample, AvgP) %>% 
  select(-gene_name) %>% 
  ggcorr(nbreaks = 10, label = F, hjust = 1, size = 3, layout.exp = 3) +
  ggtitle("correlation between mutants based on partner interaction probabilities")



saint_and_interfaces %>% 
  select(sample, gene_name, AvgP) %>% 
  spread(gene_name, AvgP) %>% 
  select(-sample) %>% 
  ggcorr(nbreaks = 10, label = T, hjust = 1, size = 3, layout.exp = 3)+
  ggtitle("correlations between core partners based on probabilites of interactions")

  


saint %>% 
  filter(gene_name %in% core_partners & BFDR < 0.05) %>% 
  select(sample, gene_name, AvgP) %>% 
  spread(gene_name, AvgP) %>% 
  select(-sample) %>% 
  ggcorr(nbreaks = 10, label = T, hjust = 1, size = 3, layout.exp = 3)+
  ggtitle("correlations between core partners based on (significant, BFDR < 0.05) probabilites of interactions")
 

saint_real <-  saint %>%
  select(Tag, sample, gene_name, AvgP, BFDR) %>% 
  filter((Tag == "C" & BFDR < 0.05) | (Tag == "N" & BFDR < 0.01)) %>%
  arrange(Tag)
saint_real_C <- saint_real %>% filter(Tag == "C")
saint_real_N <- saint_real %>% filter(Tag == "N")


saint_real %>% 
  select(Tag, sample, gene_name, AvgP, BFDR) %>% 
  arrange(Tag) %>% 
  write_tsv(path = "SAINT_intensity_cytoscape_network.txt") 

saint %>%
  filter(gene_name %in% core_partners & BFDR < 0.05) %>% 
  select(sample, gene_name, BFDR, AvgP) %>% 
  write_tsv(path = "SAINT_intensity_core_partners_cyto.txt")

preys <- saint_real %>% 
  filter(! gene_name %in% core_partners) %>% pull(gene_name) %>% unique()
saint_cyto_att <- tibble("Gene" = preys, "Type" = "Prey") %>% 
  bind_rows(., tibble("Gene" = samples_by_residue_number, "Type" = "Bait")) %>% 
  bind_rows(., tibble("Gene" = core_partners, "Type" = "Core")) %>% 
  write_tsv(path = "SAINT_intensity_cyto_att.txt")

