# Tidy Data 
# Combine plate and replicate data
library(dplyr)
library(reshape2)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(mgcv)
library(viridis)
library(pheatmap)
library(RColorBrewer)

plate1_data = read.delim('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate1-run1.txt') 
LFC_plate1_data = plate1_data %>% 
  rename('pDNA' = grep('pDNA', colnames(plate1_data), value = TRUE)) %>%
  melt(id.vars = c('Construct.Barcode', 'Construct.IDs', 'pDNA')) %>%
  mutate(LFC = value - pDNA) %>%
  select(-value, -pDNA) %>%
  spread(variable, LFC)

plate2_data = read.delim('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate2-run1.txt')
LFC_plate2_data = plate2_data %>% 
  rename('pDNA' = grep('pDNA', colnames(plate2_data), value = TRUE)) %>%
  melt(id.vars = c('Construct.Barcode', 'Construct.IDs', 'pDNA')) %>%
  mutate(LFC = value - pDNA) %>%
  select(-value, -pDNA) %>%
  spread(variable, LFC)

all_spread = inner_join(LFC_plate1_data, LFC_plate2_data, by = c('Construct.Barcode', 'Construct.IDs'))

pheatmap(all_spread %>% select(-Construct.Barcode, -Construct.IDs) %>% cor())

all_melted = all_spread %>%
  melt(id.vars = c('Construct.Barcode', 'Construct.IDs')) %>%
  separate(variable, c('Cell', 'Rep', 'Passage', 'Vector'), sep = '\\.') %>%
  group_by(Construct.Barcode, Construct.IDs, Cell, Passage, Vector) %>%
  summarise(Avg.LFC = mean(value)) %>%
  separate(Construct.IDs, c('gene1', 'gene2', 'guide1', 'guide2'), sep = ';|:')

Base_LFC = all_melted %>% filter(gene1 == 'non-targeting' |
                                     gene2 == 'non-targeting') %>%
  mutate(non_control = ifelse(gene1 == 'non-targeting', guide2, guide1)) %>%
  group_by(non_control, Cell, Passage, Vector) %>%
  summarise(base_LFC = mean(Avg.LFC), n_base = n()) 

sum_LFC = all_melted %>%
  filter(gene1 != 'non-targeting' & gene2 != 'non-targeting') %>%
  left_join(Base_LFC, by = c('guide1' = 'non_control', 'Cell', 'Passage', 'Vector')) %>%
  left_join(Base_LFC, by = c('guide2' = 'non_control', 'Cell', 'Passage', 'Vector'), 
            suffix = c('.1','.2')) %>%
  mutate(sum_LFC = base_LFC.1 + base_LFC.2) %>%
  group_by(Cell, Passage) %>%
  mutate(residual = gam(Avg.LFC ~ s(sum_LFC))$residual,
         assay = paste(Cell, Passage, sep = '_')) %>%
  ungroup()

ggplot(sum_LFC) +
  aes(x = sum_LFC, y = Avg.LFC, color = residual) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  theme(aspect.ratio = 1) +
  facet_wrap(c('Cell', 'Passage')) +
  scale_color_viridis(option = "B")

# Generate p-values for gene_combo residuals from null distribution
gene_df = sum_LFC %>%
  select(gene1, gene2) %>%
  distinct() %>%
  as.matrix() %>%
  t() %>%
  apply(2, sort) %>%
  t() %>%
  as.data.frame() %>%
  distinct() 

n_combos = nrow(gene_df)
genes = unique(sum_LFC$gene1)
assays = unique(sum_LFC$assay)
assay_pvalues = list()

for (i in 1:length(assays)) {
  curr_assay = assays[i]
  print(curr_assay)
  assay_df = sum_LFC %>%
    filter(assay == curr_assay)
  assay_residuals = assay_df$residual
  pvalues = vector(mode = 'numeric', length = n_combos)
  mean_residuals = pvalues
  mean_sum_LFCs = pvalues
  for (j in 1:n_combos) {
    if(j %% 100 == 0) {
      print(j/n_combos)
    }
    gene_V1 = gene_df$V1[j]
    gene_V2 = gene_df$V2[j]
    combo_df = assay_df %>%
      filter((gene1 == gene_V1 & gene2 == gene_V2) |
               (gene2 == gene_V1 & gene1 == gene_V2)) 
    combo_residuals = combo_df %>%
      select(residual) %>%
      unlist() 
    combo_sum_LFC = combo_df %>%
      select(sum_LFC) %>%
      unlist() 
    mean_residuals[j] = mean(combo_residuals)
    mean_sum_LFCs[j] = mean(combo_sum_LFC)
    pvalues[j] = t.test(combo_residuals, assay_residuals[-match(combo_residuals, assay_residuals)])$p.value
  }
  assay_pvalue = bind_cols(gene_df, data.frame(p_value = pvalues, 
                                               mean_residual = mean_residuals,
                                               mean_sum_LFC = mean_sum_LFCs,
                                               assay = curr_assay))
  assay_pvalues[[i]] = assay_pvalue
}
p_cutoff = 0.001

assay_pvalues_df = bind_rows(assay_pvalues) %>%
  rename(gene1 = V1, gene2 = V2) %>%
  mutate(directional_p = sign(mean_residual)*-log10(p_value),
         combo = paste(gene1, gene2, sep = '_'))

spread_df = assay_pvalues_df %>%
  group_by(combo) %>%
  filter(min(p_value) < p_cutoff) %>%
  ungroup() %>%
  select(-p_value, -mean_residual) %>%
  spread(assay, directional_p) 

spread_mat = spread_df %>%
  select(A549_P1, A549_P6, PANC1_P2, PANC1_P6, PANC1_P7) %>%
  as.matrix()

row.names(spread_mat) = spread_df$combo

pheatmap(spread_mat, fontsize = 12,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         breaks = c(seq(min(assay_pvalues_df$directional_p),0, length.out = 50),
                    seq(max(assay_pvalues_df$directional_p)/100,max(assay_pvalues_df$directional_p), length.out = 50)),
         cutree_cols =2,
         cutree_rows = 6,
         border_color = 'black',
         width = 10, height = 3,
         cluster_rows = TRUE)

ggplot(sum_LFC) +
  aes(x = sum_LFC, y = Avg.LFC, color = residual) +
  geom_point(aes(alpha = )) +
  theme(aspect.ratio = 1) +
  facet_wrap(c('Cell', 'Passage')) +
  scale_color_viridis(option = "B")

gene_j = 'AP1M1'
gene_i = 'AP1M2'

interesting_LFC = sum_LFC %>%
  mutate(interesting = ((gene1 == gene_j) & (gene2 == gene_i)) |
           ((gene2 == gene_j) & (gene1 == gene_i))) %>%
  arrange(interesting)

ggplot(interesting_LFC) +
  aes(x = sum_LFC, y = Avg.LFC) +
  geom_point(aes(color = interesting, shape = interesting, size = interesting, alpha = interesting)) +
  scale_size_manual(values = c(1,2))+
  scale_alpha_manual(values = c(0.1, 0.75)) +
  scale_color_brewer(palette = 'Paired') +
  theme(aspect.ratio = 1) +
  theme_classic() +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  facet_wrap('assay')

ggplot(assay_pvalues_df) +
  aes(x = mean_sum_LFC, y = directional_p) +
  geom_point() +
  theme(aspect.ratio = 1) +
  facet_wrap('assay')
