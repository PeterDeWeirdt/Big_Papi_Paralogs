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

get_base_LFC <- function(melted, FUN) {
  Base_LFC = melted %>% 
    filter(gene1 == 'non-targeting' |
             gene2 == 'non-targeting') %>%
    mutate(non_control = ifelse(gene1 == 'non-targeting', guide2, guide1)) %>%
    group_by(non_control, Cell, Passage, Vector) %>%
    summarise(base_LFC = FUN(Avg.LFC), n_base = n()) 
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

summarise_combos <- function(guide_residuals, FUN) {
  gene_residuals = guide_residuals %>%
    mutate(temp_gene1 = ifelse(gene1 < gene2, gene1, gene2),
           temp_gene2 = ifelse(gene1 < gene2, gene2, gene1)) %>%
    group_by(temp_gene1, temp_gene2, assay) %>%
    summarise(Cell = first(Cell),
              Passage = first(Passage),
              Combo.LFC = FUN(Avg.LFC),
              Gene.LFC.1 = FUN(base_LFC.1), 
              Gene.LFC.2 = FUN(base_LFC.2),
              Combo.base.LFC = FUN(sum_LFC),
              Combo.residual = FUN(residual)) %>%
    rename(gene1 = temp_gene1, gene2 = temp_gene2) %>%
    ungroup() 
  return(gene_residuals)
}

generate_pvalues <- function(calculated.values, calc.cols, original.values, orig.cols, my.stat, n.samps,
                             nboot = 1e4, group.col = 'assay') {
  # calculated.values: df with (1) statistical values and (2) grouping variable. 
  # original.values: similaar to calculated.valuess, but not summarised
  calculated.values = calculated.values %>% select(calc.cols)
  original.values = original.values %>% select(orig.cols)
  assays = unique(calculated.values[,2])
  all_pvals = c()
  for (curr_assay in assays) {
    print(curr_assay)
    values = original.values %>% filter(!!as.name(group.col) == curr_assay) %>% .[,1] %>% unlist()
    stats = vector(mode = 'numeric', length = nboot)
    for (i in 1:nboot) {
      if(i %% 10000 == 0) {print(paste('bootstrap progress:', as.character(i/nboot)))}
      stats[i] = my.stat(sample(values, n.samps))
    }
    
    calc.stats = calculated.values %>% filter(!!as.name(group.col) == curr_assay) %>% .[,1] %>% unlist()
    pvals = vector(mode = 'numeric', length = length(calc.stats))
    for (i in 1:length(pvals)) {
      pvals[i] = min(sum(calc.stats[i] > stats), sum(calc.stats[i] < stats))/
        length(stats)
    }
    all_pvals = c(pvals, all_pvals)
  }
  return(all_pvals)
}

citations_helper <- function(row, i, denom) {
  # Query the biogrid database for citations
  # Helpful querrying terms can be found - https://wiki.thebiogrid.org/doku.php/biogridrest
  if (i %% 10 == 0) {print(i/denom)}
  citations = read.delim(url(paste('https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=',
                                      row[[1]],'|',
                                      row[[2]],
                                      '&taxId=9606&accesskey=ef396dc62ddfc51b05547bf9dd6feb32&includeInteractors=false&includeHeader=true&format=tab2&selfcitationsExcluded=true', 
                                      sep = '')))
  n = nrow(citations)
}

get_interactors <- function(gene_combos) {
  n_citations = map_int(1:nrow(gene_combos), 
                                  function(i) citations_helper(gene_combos[i,], i, nrow(gene_combos)))
  gene_combos['citations'] = n_citations
  return(gene_combos)
}

my_ks <- function(value, condition) {
  p = ks.test(value[condition > 0], value[condition == 0], alternative = 'two.sided')$p.value
  return(paste('KS p:', as.character(signif(p, 3))))
}

plot_rug <- function(residual_biogrid) {
  residual_biogrid = ungroup(residual_biogrid)
  p = ggplot(residual_biogrid) +
    aes(x = Combo.residual, color = citations > 0) +
    stat_ecdf() + 
    facet_wrap('assay') + 
    geom_rug(sides = 'b') +
    scale_color_brewer(palette = 'Paired') +
    base_theme() +
    geom_text(aes(label = label,  x = m.x, y = m.y, hjust = m.hjust, vjust = m.vjust), 
              data = residual_biogrid %>%
                group_by(assay) %>%
                summarise(label = my_ks(Combo.residual, citations)) %>%
                mutate(m.x = -Inf, m.y = Inf, m.hjust = -0.1, m.vjust = 1.5),
              inherit.aes = FALSE)
  return(p)
}

compare_scores <- function(residuals_x, residuals_y, xlabel, ylabel, p.cut, label_p) {
  joined_residuals = inner_join(residuals_x, residuals_y, 
                                by = c('gene1', 'gene2', 'assay'), 
                                suffix = c('.x', '.y')) %>%
    mutate(biogrid = citations.x > 0)
  subset_data = joined_residuals %>% 
    group_by(assay) %>%
    summarise(y_min = max(Combo.residual.y[p.value.y < p.cut & Combo.residual.y < 0]), 
              x_min = max(Combo.residual.x[p.value.x < p.cut & Combo.residual.x < 0]), 
              y_max = min(Combo.residual.y[p.value.y < p.cut & Combo.residual.y > 0]), 
              x_max = min(Combo.residual.x[p.value.x < p.cut & Combo.residual.x > 0]))
    
  
  p = ggplot(joined_residuals) +
    aes(x = Combo.residual.x, y = Combo.residual.y) +
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
    ggtitle('Gene residuals') +
    geom_text_repel(data = joined_residuals %>% filter(gene1 == label_p[1], gene2 == label_p[2], 
                                                       assay == label_p[3]), 
                     aes(label = paste(gene1, gene2, sep = '_')), point.padding = 0.5, nudge_y = 0.1,  
                    min.segment.length = 0, size = 3) +
    stat_cor() 
    
  return(p)
}

visualize_combo <- function(melted, my.gene1, my.gene2, assay, nt = 'non-targeting') {
  # be sure gene1 and gene2 are alphabetical
  melted_combo = melted %>% 
    mutate(temp_gene1 = ifelse(gene1 == nt, gene2, ifelse(gene2 == nt, gene1, ifelse(gene1 < gene2, gene1, gene2))),
           temp_gene2 = ifelse(gene2 == nt | gene1 == nt, nt, ifelse(gene1 < gene2, gene2, gene1))) %>%
    select(-gene1, -gene2) %>%
    rename(gene1 = temp_gene1, 
           gene2 = temp_gene2) %>%
    filter(paste(Cell, Passage, sep = '_') == assay &
             ((gene1 == my.gene1 & gene2 == nt) | 
                (gene1 == my.gene2 & gene2 == nt) |
                (gene1 == my.gene1 & gene2 == my.gene2))) %>%
    mutate(combo = paste(gene1, gene2, sep='_'))
  p=ggplot(melted_combo) +
    aes(x = combo, y = Avg.LFC, color = combo) +
    base_theme() +
    geom_boxplot(alpha = 0) +
    geom_point(position = position_dodge(width = 0.7)) +
    scale_color_brewer(palette = 'Set2') +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") +
    theme(legend.position = '', axis.text.x = element_text(angle = 45,hjust = 1)) +
    ggtitle(assay)
  
  return(p)
}

create_heatmap <- function(combo_residuals,p_cut) {
  signif_combos = combo_residuals %>% 
    mutate(combo = paste(gene1, gene2, sep = '_')) %>%
    group_by(assay) %>%
    group_by(combo) %>%
    filter(min(p.value) < p_cut)
  
  combo_order_df = signif_combos %>% 
    group_by(combo) %>%
    summarise(max_abs_resid = Combo.residual[which.max(abs(Combo.residual))]) %>%
    arrange(-max_abs_resid)
  
  assay_order_df = signif_combos %>%
    group_by(assay) %>%
    summarise(P = first(Passage)) %>%
    arrange(P)
  
  signif_combos['combo'] = factor(signif_combos$combo, levels = order_df$combo)
  signif_combos['assay'] = factor(signif_combos$assay, 
                                  levels = assay_order_df$assay)
  
  ggplot(signif_combos) +
    aes(fill = Combo.residual, x = assay, y = combo) +
    geom_tile() +
    base_theme() +
    scale_fill_viridis(option = 'C') +
    theme(axis.text.x = element_text(angle = 90))
 }


