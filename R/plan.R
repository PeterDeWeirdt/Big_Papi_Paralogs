# The workflow plan data frame outlines what you are going to do.
plan <- drake_plan(
  plate1_data = read.delim(file_in('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate1-run1.txt')),
  plate2_data = read.delim(file_in('./inputs/lognorm-JD_GPP568_Raum_20190115_sgRNA_Plate2-run1.txt')),
  all_spread = clean_inputs(list(plate1_data, plate2_data)),
  rep_heatmap = pheatmap(all_spread %>% select(-Construct.Barcode, -Construct.IDs) %>% cor(), 
                         main = 'Correlation between guide LFCs'),
  all_melted = melt_spread(all_spread),
  LFC_dist = plot_lfc(all_melted),
  base_LFC = get_base_LFC(all_melted),
  sum_LFC = model_expectation(all_melted, base_LFC),
  residuals_plot = plot_residuals(sum_LFC),
  combo_residuals = summarise_combos(sum_LFC),
  biogrid_interactions = get_interactors(combo_residuals %>% select(gene1, gene2) %>% distinct()),
  residual_biogrid = inner_join(combo_residuals, biogrid_interactions, by = c('gene1', 'gene2')),
  pr = get_auc(residual_biogrid, 'mean.residual'),
  pr_plot = plot_pr(pr), 
  
  # Now do the same steps, but this time calling outliers
  cutoff = 1e-2,
  outliers = calc_outliers(all_melted),
  plot_outlier_base_dist = plot_out_base(outliers, base_LFC, cutoff = cutoff), 
  melted_outliers = join_melted_outliers(all_melted, outliers), 
  
  # method 1: just call non-targetting outliers
  nt_filtered_outliers = melted_outliers %>% group_by(Cell, Passage) %>%
    filter(!((gene1 == 'non-targeting' & (p.value.guide1 < cutoff)) |
               (gene2 == 'non-targeting' & (p.value.guide2 < cutoff)))), 
  outlier_base_LFC = get_base_LFC(nt_filtered_outliers),
  nt_sum_LFC = model_expectation(nt_filtered_outliers, outlier_base_LFC),
  nt_residuals_plot = plot_residuals(nt_sum_LFC),
  nt_combo_residuals = summarise_combos(nt_sum_LFC),
  nt_residual_biogrid = inner_join(nt_combo_residuals, biogrid_interactions, by = c('gene1', 'gene2')),
  nt_pr = get_auc(nt_residual_biogrid, 'mean.residual'),
  
  # method 2 call all outliers
  
  all_called_outliers = nt_filtered_outliers %>% 
    mutate(outlier = p.value.guide1 < cutoff | p.value.guide2 < cutoff),
  all_sum_LFC = model_expectation(all_called_outliers %>% filter(outlier), outlier_base_LFC),
  all_residuals_plot = plot_residuals(all_sum_LFC),
  all_combo_residuals = summarise_combos(all_sum_LFC),
  all_residual_biogrid = inner_join(all_combo_residuals, biogrid_interactions, by = c('gene1', 'gene2')),
  all_pr = get_auc(all_residual_biogrid, 'mean.residual'),

  # method 3: call outliers that are non-lethal and significant
  
  called_outliers = call_outliers(nt_filtered_outliers, outlier_base_LFC, cutoff),
  outlier_LFC_dist = plot_lfc(called_outliers, 'outlier'),
  outlier_sum_LFC = model_expectation(called_outliers %>% 
                                        filter(outlier != TRUE) %>%
                                        select(-p.value.guide1, -p.value.guide2,
                                               -outlier, -base_LFC.1, -base_LFC.2), outlier_base_LFC),
  outlier_residuals_plot = plot_residuals(outlier_sum_LFC),
  outlier_combo_residuals = summarise_combos(outlier_sum_LFC),
  outlier_residual_biogrid = inner_join(outlier_combo_residuals, biogrid_interactions, by = c('gene1', 'gene2')),
  outlier_pr = get_auc(outlier_residual_biogrid, 'mean.residual'),
  outlier_pr_plot = plot_pr(outlier_pr),
  
  # Compare performance of interesting scores
  joined_pr = pr %>% 
    mutate(filter = 'no filtering') %>%
    bind_rows(outlier_pr %>% mutate(filter = 'non-lethal outliers'), 
              nt_pr %>% mutate(filter = 'NT'), 
              all_pr %>% mutate(filter = 'all outliers')) %>%
    mutate(filter = factor(filter, levels = c('no filtering','NT', 'non-lethal outliers','all outliers'))),
  joined_pr_plot = plot_pr(joined_pr, 'filter'),
  
  # Heatmaps 
  nf_heat = plot_hit_heatmaps(residual_biogrid, 'No filtering hits'),
  filt_heat = plot_hit_heatmaps(all_residual_biogrid, 'Filter all outliers'),
  
  # Comparing scores w/ scatter
  interaction_mat = join_interactions(residual_biogrid, all_residual_biogrid, '.no_filter', '.all_filter'),
  interaction_output = paste(Sys.Date(), 'filtering_scheme_interactions.csv', sep = '_'),
  write.csv(interaction_mat, file_out("2019-02-26_filtering_scheme_interactions.csv")),
  
  scatter_comp = compare_scores(residual_biogrid, all_residual_biogrid, 'no filtering', 'filter all outliers'),
  
  # Send plots to rmarkdown 
  knitr_output = paste(Sys.Date(), "outlier_report.html", sep = '_'),
  report = rmarkdown::render(
    knitr_in("outlier_report.Rmd"),
    output_file = file_out("2019-02-26_outlier_report.html"),
    quiet = TRUE
  )
)
