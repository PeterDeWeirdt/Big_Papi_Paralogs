# The workflow plan data frame outlines what you are going to do.
plan <- drake_plan(
  plate1_data = read.delim(file_in('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate1-run1.txt')),
  plate2_data = read.delim(file_in('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate2-run1.txt')),
  all_spread = clean_inputs(list(plate1_data, plate2_data)),
  rep_heatmap = pheatmap(all_spread %>% select(-Construct.Barcode, -Construct.IDs) %>% cor(), 
                         main = 'Correlation between guide LFCs'),
  all_melted = melt_spread(all_spread),
  biogrid_interactions = get_interactors(all_melted %>% select(gene1, gene2) %>% distinct()),
  LFC_dist = plot_lfc(all_melted),
  
  
  # mean ressiduals
  base_LFC_mean = get_base_LFC(all_melted, mean),
  sum_LFC_mean = model_expectation(all_melted, base_LFC_mean),
  residuals_plot_mean = plot_residuals(sum_LFC_mean),
  combo_residuals_mean = summarise_combos(sum_LFC_mean, mean),
  pvalues_mean = generate_pvalues(combo_residuals_mean, c('Combo.residual', 'assay'), 
                                  sum_LFC_mean, c('residual', 'assay'), 
                                  mean, 4, 1e4),
  combo_pvalues_mean = combo_residuals_mean %>% mutate(p.value = pvalues_mean),
  residual_biogrid_mean = inner_join(combo_pvalues_mean, biogrid_interactions, by = c('gene1', 'gene2')),
  ks_plot_mean = plot_rug(residual_biogrid_mean),
  
  # median residuals
  base_LFC_median = get_base_LFC(all_melted, median),
  sum_LFC_median = model_expectation(all_melted, base_LFC_median),
  residuals_plot_median = plot_residuals(sum_LFC_median),
  combo_residuals_median = summarise_combos(sum_LFC_median, median),
  pvalues_median = generate_pvalues(combo_residuals_median, c('Combo.residual', 'assay'), 
                                  sum_LFC_median, c('residual', 'assay'), 
                                  median, 4, 1e4),
  combo_pvalues_median = combo_residuals_median %>% mutate(p.value = pvalues_median),
  residual_biogrid_median = inner_join(combo_pvalues_median, biogrid_interactions, by = c('gene1', 'gene2')),
  ks_plot_median = plot_rug(residual_biogrid_median),
  
  # compare mean to median
  comparison = compare_scores(residual_biogrid_mean, residual_biogrid_median, 'mean', 'median', 5e-2,
                              c('SEC23A', 'SEC23B', 'PANC1_P7')),
  
  # Move forward with mean
  SEC_p = visualize_combo(all_melted, 'SEC23A', 'SEC23B', 'PANC1_P7'),
  MAPK_p = visualize_combo(all_melted, 'MAPK1', 'MAPK3', 'PANC1_P7'),
  
  # Create heatmap from mean data
  heatmap = create_heatmap(residual_biogrid_mean, 0.05)
  
)
