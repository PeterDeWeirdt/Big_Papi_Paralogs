# Load our packages in the usual way.
# This example requires packages forcats, readxl, and rmarkdown,
# but you do not need to load them here.
library(drake)
require(dplyr)
require(ggplot2)
library(tidyr)
library(purrr)
library(ClassDiscovery)
library(mgcv)
library(splines)
library(PRROC)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
pkgconfig::set_config("drake::strings_in_dots" = "literals") # New file API


base_theme <- function() {
  #Provide base elements for ggplot
  return(theme_classic() +
    theme(text = element_text(size = 12),
          aspect.ratio = 1)
  )
}

clean_inputs <- function(data_list) {
  clean = list()
  for(i in 1:length(data_list)) {
    data = data_list[[i]]
    clean[[i]] = data %>%
      rename('pDNA' = grep('pDNA', colnames(data), value = TRUE)) %>%
      melt(id.vars = c('Construct.Barcode', 'Construct.IDs', 'pDNA')) %>%
      mutate(LFC = value - pDNA) %>%
      select(-value, -pDNA) %>%
      spread(variable, LFC)
  }
  all_spread = plyr::join_all(clean, by = c('Construct.Barcode', 'Construct.IDs'), type = 'inner')
  return(all_spread)
}

melt_spread <- function(data) {
  all_melted = data %>%
    melt(id.vars = c('Construct.Barcode', 'Construct.IDs')) %>%
    separate(variable, c('Cell', 'Rep', 'Passage', 'Vector'), sep = '\\.') %>%
    group_by(Construct.Barcode, Construct.IDs, Cell, Passage, Vector) %>%
    summarise(Avg.LFC = mean(value)) %>%
    separate(Construct.IDs, c('gene1', 'gene2', 'guide1', 'guide2'), sep = ';|:') %>%
    ungroup()
  # Put genes in alphabetical order
  return(all_melted)
}

plot_lfc <- function(melted, color_var = 'non-targeting') {
  p = ggplot(melted %>%
               mutate(curr_color = ifelse(color_var == 'non-targeting', 
                                            ((color_var == gene1) | (color_var == gene2)),
                                            !!as.name(color_var)))) +
    aes(x = Avg.LFC, color = outlier) +
    geom_density() +
    facet_wrap(c('Cell', 'Passage')) +
    ggtitle(paste('LFC distirubtion compared with', color_var)) +
    base_theme() +
    labs(color = color_var)
  return(p)
}

get_base_LFC <- function(melted) {
  Base_LFC = melted %>% 
    filter(gene1 == 'non-targeting' |
             gene2 == 'non-targeting') %>%
    mutate(non_control = ifelse(gene1 == 'non-targeting', guide2, guide1)) %>%
    group_by(non_control, Cell, Passage, Vector) %>%
    summarise(base_LFC = mean(Avg.LFC), n_base = n()) 
  return(Base_LFC)
}

model_expectation <- function(melted, Base_LFC) {
  # model the activity of combos as the sum of their base LFC's
  sum_LFC = melted %>%
    filter(gene1 != 'non-targeting' & gene2 != 'non-targeting') %>%
    left_join(Base_LFC, by = c('guide1' = 'non_control', 'Cell', 'Passage', 'Vector')) %>%
    left_join(Base_LFC, by = c('guide2' = 'non_control', 'Cell', 'Passage', 'Vector'), 
              suffix = c('.1','.2')) %>%
    mutate(sum_LFC = base_LFC.1 + base_LFC.2) %>%
    group_by(Cell, Passage) %>%
    mutate(residual = gam(Avg.LFC ~ ns(sum_LFC, df = 3))$residual,
           assay = paste(Cell, Passage, sep = '_')) %>%
    ungroup()
  return(sum_LFC)
}

plot_residuals <- function(modeled_data) {
  ggplot(modeled_data) +
    aes(x = sum_LFC, y = Avg.LFC, color = residual) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam", formula = y ~ ns(x, df = 3)) +
    geom_smooth(method = 'lm', formula = y ~ offset(x), color = 'black', size = 0.5) +
    facet_wrap(c('Cell', 'Passage')) +
    scale_color_viridis(option = "B") +
    base_theme()
}

summarise_combos <- function(guide_residuals) {
  gene_residuals = guide_residuals %>%
    mutate(temp_gene1 = ifelse(gene1 < gene2, gene1, gene2),
           temp_gene2 = ifelse(gene1 < gene2, gene2, gene1)) %>%
    group_by(temp_gene1, temp_gene2, assay) %>%
    summarise(Cell = first(Cell),
              Passage = first(Passage),
              Avg.LFC = mean(Avg.LFC),
              base.Avg.LFC.1 = mean(base_LFC.1), 
              base.Avg.LFC.2 = mean(base_LFC.2),
              sum.base.LFC = mean(sum_LFC),
              mean.residual = mean(residual)) %>%
    rename(gene1 = temp_gene1, gene2 = temp_gene2) 
  return(gene_residuals)
}

interactions_helper <- function(row, i, denom) {
  # Query the biogrid database for interactions
  # Helpful querrying terms can be found - https://wiki.thebiogrid.org/doku.php/biogridrest
  if (i %% 10 == 0) {print(i/denom)}
  interactions = read.delim(url(paste('https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=',
                                      row[[1]],'|',
                                      row[[2]],
                                      '&taxId=9606&accesskey=ef396dc62ddfc51b05547bf9dd6feb32&includeInteractors=false&includeHeader=true&format=tab2&selfInteractionsExcluded=true', 
                                      sep = '')))
  n = nrow(interactions)
}

get_interactors <- function(gene_combos) {
  n_interactions = map_int(1:nrow(gene_combos), 
                                  function(i) interactions_helper(gene_combos[i,], i, nrow(gene_combos)))
  gene_combos['interactions'] = n_interactions
  return(gene_combos)
}

auc_helper <- function(ranks, interactions) {
  positive_ranks = ranks[interactions != 0]
  negative_ranks = ranks[interactions == 0]
  pr = pr.curve(positive_ranks, negative_ranks,
                dg.compute = FALSE)
  return(pr$auc.integral)
}

get_auc <- function(combo_interactions, score_column) {
  assay_aucs = combo_interactions %>%
    group_by(assay) %>%
    summarise(auc = auc_helper(abs(!!as.name(score_column)), interactions), 
              random = sum(interactions != 0)/n(),
              auc.ratio = auc/random)
}

plot_pr <- function(pr_data, fill = 'assay') {
  p = ggplot(pr_data) +
    aes(x = assay, y = auc.ratio, fill = !!as.name(fill)) +
    geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
    scale_fill_brewer(palette = 'Set2') +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    base_theme() +
    theme(axis.text.x = element_text(angle = 90)) 
    
  return(p)
}

outlier_helper <- function(gene_data, dominant.guide, secondary.guide) {
  spread_guides = gene_data %>% 
    select(!!as.name(dominant.guide), !!as.name(secondary.guide), Avg.LFC) %>%
    spread(!!as.name(secondary.guide), Avg.LFC) %>%
    select_if(~ !any(is.na(.))) # If some guides not tested in all conditions
  spread_mat = as.matrix(spread_guides[,-1])
  row.names(spread_mat) = spread_guides[[dominant.guide]]
  spca <- SamplePCA(t(spread_mat))
  mQC = mahalanobisQC(spca, 2)
  mQC['guide'] = spread_guides[[dominant.guide]]
  return(mQC)
}

calc_outliers <- function(melted_data) {
  melted_data = melted_data %>% ungroup()
  assay_list = list()
  for (cell in (unique(melted_data$Cell))) {
    cell_data = melted_data %>%
      filter(Cell == cell)
    for (passage in unique(cell_data$Passage)) {
      assay_data = cell_data %>%
        filter(Passage == passage) 
      print(paste(cell, passage))
      gene_list = list()
      for(curr_gene in unique(assay_data$gene1)) {
        gene1_outliers = outlier_helper(assay_data %>% 
                                          filter(gene1 == curr_gene),
                                        'guide1', 'guide2') %>%
          mutate(dominant.gene = 'gene1',
                 curr.gene = curr_gene)
        gene2_outliers = outlier_helper(assay_data %>% 
                                          filter(gene2 == curr_gene), 
                                        'guide2', 'guide1') %>%
          mutate(dominant.gene = 'gene2',
                 curr.gene = curr_gene)
        gene_list[[curr_gene]] = bind_rows(gene1_outliers, gene2_outliers)
      }
      assay_list[[paste(cell, passage, sep = '_')]] = bind_rows(gene_list) %>%
        mutate(Cell = cell, Passage = passage)
    }
  }
  return(bind_rows(assay_list))
}

join_melted_outliers <-function(all_melted, outliers) {
  joined = left_join(all_melted, outliers %>% 
              select(p.value, guide, Cell, Passage) %>%
              rename(p.value.guide1 = p.value), 
            by = c('Cell', 'Passage', 'guide1' = 'guide')) %>%
    left_join(outliers %>% 
                select(p.value, guide, Cell, Passage) %>%
                rename(p.value.guide2 = p.value), 
              by = c('Cell', 'Passage', 'guide2' = 'guide'))
  return(joined)
}

call_outliers <- function(melted_outliers, Base_LFC, cutoff) {
  based_outliers = melted_outliers %>%
    left_join(Base_LFC, by = c('guide1' = 'non_control', 'Cell', 'Passage', 'Vector')) %>%
    left_join(Base_LFC, by = c('guide2' = 'non_control', 'Cell', 'Passage', 'Vector'), 
              suffix = c('.1','.2')) %>%
    group_by(Cell, Passage, gene1) %>%
    mutate(gene1.buff = base_LFC.1 > mean(base_LFC.1), 
           gene1.leth = base_LFC.1 < mean(base_LFC.1)) %>%
    group_by(Cell, Passage, gene2) %>%
    mutate(gene2.buff = base_LFC.2 > mean(base_LFC.2), 
           gene2.leth = base_LFC.2 < mean(base_LFC.2)) %>%
    group_by(Cell, Passage) %>%
    mutate(buff.outlier = (gene1.buff & (p.value.guide1 < cutoff)) | 
             (gene2.buff & (p.value.guide2 < cutoff)), 
           leth.outlier = (gene1.leth & (p.value.guide1 < cutoff)) | 
             (gene2.leth & (p.value.guide2 < cutoff))) %>%
    group_by(Cell, Passage, gene1) %>%
    mutate(gene1.nonleth.outlier = (p.value.guide1 > cutoff) & (max(gene1.leth & (p.value.guide1 < cutoff)) == 1)) %>%
    group_by(Cell, Passage, gene2) %>%
    mutate(gene2.nonleth.outlier = (p.value.guide2 > cutoff) & (max(gene2.leth & (p.value.guide2 < cutoff)) == 1)) %>%
    mutate(outlier = buff.outlier | gene1.nonleth.outlier | gene2.nonleth.outlier) %>%
    ungroup()
  return(based_outliers)
}

plot_out_base <- function(outliers, base, cutoff) {
  joined = inner_join(outliers, base, by = c('Cell', 'Passage', 'guide' = 'non_control')) %>%
    mutate(outlier = p.value < cutoff)
  ggplot(joined) +
    aes(x = base_LFC, fill = outlier) +
    geom_density(alpha = 0.5) +
    facet_wrap(c('Cell', 'Passage')) +
    ggtitle(paste('Base LFC distirubtion compared with outliers')) +
    base_theme() 
}

plot_hit_heatmaps <- function(residual_df, title) {
  annotation_df = residual_df %>% 
    ungroup() %>%
    mutate(interaction = paste(gene1, gene2, sep = '_'), 
           biogrid = as.character(interactions > 0)) %>%
    select(interaction, biogrid) %>%
    distinct() %>%
    as.data.frame()
  interactions = annotation_df$interaction
  annotation_df = annotation_df %>% select(biogrid)
  row.names(annotation_df) = interactions
  
  spread_df = residual_df %>%
    mutate(interaction = paste(gene1, gene2, sep = '_')) %>%
    group_by(assay) %>%
    mutate(hit =  abs(mean.residual - mean(mean.residual)) > 2*sd(mean.residual)) %>%
    group_by(interaction) %>%
    filter(max(hit) == 1) %>%
    ungroup() %>%
    select(interaction, assay, mean.residual) %>%
    spread(assay, mean.residual) 
  
  spread_mat = spread_df %>%
    select(A549_P1, A549_P6, PANC1_P2, PANC1_P6, PANC1_P7) %>%
    as.matrix()
  
  row.names(spread_mat) = spread_df$interaction
  
  p = pheatmap(spread_mat,
           color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "RdBu")))(100),
           breaks = c(seq(min(residual_df$mean.residual),0, length.out = 50),
                      seq(max(residual_df$mean.residual)/100,max(residual_df$mean.residual), length.out = 50)),
           cutree_cols =2,
           cutree_rows = 6,
           border_color = 'black',
           width = 10, height = 3,
           cluster_rows = TRUE,
           annotation_row = annotation_df, 
           main = title, 
           annotation_colors = 
             list(biogrid = c('TRUE' = '#a6cee3',
                              'FALSE' = 'white')))
  return(p)
}

compare_scores <- function(residuals_x, residuals_y, xlabel, ylabel) {
  joined_residuals = inner_join(residuals_x, residuals_y, 
                                by = c('gene1', 'gene2', 'assay'), 
                                suffix = c('.x', '.y')) %>%
    mutate(biogrid = interactions.x > 0)
  subset_data = joined_residuals %>% 
    group_by(assay) %>%
    summarise(y_min = mean(mean.residual.y) - 2*sd(mean.residual.y), 
              x_min = mean(mean.residual.x) - 2*sd(mean.residual.x), 
              y_max = mean(mean.residual.y) + 2*sd(mean.residual.y), 
              x_max = mean(mean.residual.x) + 2*sd(mean.residual.x))
    
  
  p = ggplot(joined_residuals) +
    aes(x = mean.residual.x, y = mean.residual.y) +
    facet_wrap('assay') +
    geom_point(aes(color = biogrid)) +
    scale_color_brewer(palette = 'Paired') +
    base_theme() +
    xlab(xlabel) +
    ylab(ylabel) +
    geom_hline(data = subset_data, aes(yintercept = y_min), linetype = 'dotted') +
    geom_hline(data = subset_data, aes(yintercept = y_max), linetype = 'dotted') +
    geom_vline(data = subset_data, aes(xintercept = x_max), linetype = 'dotted') +
    geom_vline(data = subset_data, aes(xintercept = x_min), linetype = 'dotted') +
    ggtitle('Mean Residuals for Interactions') +
    stat_cor() 
    
  return(p)
}

join_interactions <- function(residuals_x, residuals_y, xlabel, ylabel) {
  joined_residuals = inner_join(residuals_x %>%
                                  group_by(assay) %>%
                                  mutate(hit = abs(mean.residual - mean(mean.residual)) > 2*sd(mean.residual)), 
                                residuals_y %>%
                                  group_by(assay) %>%
                                  mutate(hit = abs(mean.residual - mean(mean.residual)) > 2*sd(mean.residual)), 
                                by = c('gene1', 'gene2', 'assay', 'Cell', 'Passage', 'interactions'), 
                                suffix = c(xlabel, ylabel))
  return(joined_residuals)
    
}


