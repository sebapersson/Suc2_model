# This script creates the diagnostic plots for nlme. It will also output output that can be copied to 
# LaTex, to get this output run the script in rstudio (or directly in R) without piping the output. This script 
# uses the function in Process_monolix_result.R. 

source("Process_monolix_result.R")

# -----------------------------------------------------------------------------------------------------------------
# Model 1
# -----------------------------------------------------------------------------------------------------------------
# Shrinkage
path_to_result <- "../Monolix/Model1/Result"
estimated_sigma <- 0.43233
n_param <- 12
calculate_shrinkage_model(path_to_result, estimated_sigma, n_param)

# Obs vs ipred
path_obs_vs_ipred = "../../Result/Figures/Model1_nlme/Model1_obs_vs_IPRED.pdf"
print_obs_vs_ipred(path_obs_vs_ipred, path_to_result, estimated_sigma, print_figure = T)

# Residuals
path_res_time <- "../../Result/Figures/Model1_nlme/Model1_res_time.pdf"
path_res_response <- "../../Result/Figures/Model1_nlme/Model1_res_response.pdf"
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = T)
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = F, save_time = T)

# Check distribtion 
path_hist_IWRES <- "../../Result/Figures/Model1_nlme/Model1_hist_IWRES.pdf"
path_qq_IWRES <- "../../Result/Figures/Model1_nlme/Model1_qq_IWRES.pdf"
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, print_figure = T)
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, save_hist = T)

# Check individual parameters 
path_hist_ind <- "../../Result/Figures/Model1_nlme/Model1_hist_ind.pdf"
path_qq_ind <- "../../Result/Figures/Model1_nlme/Model1_qq_ind.pdf"
param_names <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc20", "Glc0")
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            n_col = 6, n_param = 12)
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            save_hist = F, n_col = 6, n_param = 12)

# Check correlation individual parameters 
path_corrplot_ind <- "../../Result/Figures/Model1_nlme/Model1_corrplot_ind.pdf"
path_pair_plot_ind <- "../../Result/Figures/Model1_nlme/Model1_pairs_ind.pdf"
plot_list_corr <- check_correlation_ind_param(path_corrplot_ind, path_pair_plot_ind, path_to_result)
pdf(path_corrplot_ind, width = 9, height = 5); corrplot.mixed(plot_list_corr[[1]]); garbage <- dev.off()
pdf(path_pair_plot_ind, width = 9, height = 5); plot_list_corr[[2]]; garbage <- dev.off()

# Check individual fits 
path_ind_fit <- "../../Result/Figures/Model1_nlme/Model1_ind_fit_"
plot_individual_fit(path_ind_fit, path_to_result, 74, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 73, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 10, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 26, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 55, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 1, print_fig = T, save_fig = T)

# Print the likelihood values 
print_likelihood_values(path_to_result)

# Print the fisher correlation matrix
path_fish_cor_mat <- "../../Result/Figures/Model1_nlme/Model1_fish_cor_mat.pdf"
plot_heat_map_fisher(path_fish_cor_mat, path_to_result, print_figure = F)

# Plot the small correlation matrix 
path_fish_cor_mat_small <- "../../Result/Figures/Model1_nlme/Model11_fish_cor_mat_small.pdf"
cor_mat <- plot_heat_map_fisher_small(path_fish_cor_mat_small, path_to_result, print_figure = T)

# Table for the population parameters 
print(create_latex_table_pop_param(path_to_result), include.rownames = F)

# Summary of the EBE:s
print_summary_stat_EBE(str_c(path_to_result, "/"))

# Coefficents of variation 
calc_coeff_of_variation(path_to_result)

# -----------------------------------------------------------------------------------------------------------------
# Model 2
# -----------------------------------------------------------------------------------------------------------------
# Shrinkage
path_to_result <- "../Monolix/Model2/Result"
estimated_sigma <- 0.4792
n_param <- 12
calculate_shrinkage_model(path_to_result, estimated_sigma, n_param)

# Obs vs ipred
path_obs_vs_ipred = "../../Result/Figures/Model2_nlme/Model2_obs_vs_IPRED.pdf"
print_obs_vs_ipred(path_obs_vs_ipred, path_to_result, estimated_sigma, print_figure = T)

# Residuals
path_res_time <- "../../Result/Figures/Model2_nlme/Model2_res_time.pdf"
path_res_response <- "../../Result/Figures/Model2_nlme/Model2_res_response.pdf"
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = T)
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = F, save_time = T)

# Check distribtion 
path_hist_IWRES <- "../../Result/Figures/Model2_nlme/Model2_hist_IWRES.pdf"
path_qq_IWRES <- "../../Result/Figures/Model2_nlme/Model2_qq_IWRES.pdf"
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, print_figure = T)
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, save_hist = T)

# Check individual parameters 
path_hist_ind <- "../../Result/Figures/Model2_nlme/Model2_hist_ind.pdf"
path_qq_ind <- "../../Result/Figures/Model2_nlme/Model2_qq_ind.pdf"
param_names <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc20", "Glc0")
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            n_col = 6, n_param = 12)
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            save_hist = F, n_col = 6, n_param = 12)

# Check correlation individual parameters 
path_corrplot_ind <- "../../Result/Figures/Model2_nlme/Model2_corrplot_ind.pdf"
path_pair_plot_ind <- "../../Result/Figures/Model2_nlme/Model2_pairs_ind.pdf"
plot_list_corr <- check_correlation_ind_param(path_corrplot_ind, path_pair_plot_ind, path_to_result)
pdf(path_corrplot_ind, width = 9, height = 5); corrplot.mixed(plot_list_corr[[1]]); garbage <- dev.off()
pdf(path_pair_plot_ind, width = 9, height = 5); plot_list_corr[[2]]; garbage <- dev.off()

# Check individual fits 
path_ind_fit <- "../../Result/Figures/Model2_nlme/Model2_ind_fit_"
plot_individual_fit(path_ind_fit, path_to_result, 74, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 73, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 10, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 26, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 55, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 1, print_fig = T, save_fig = T)

# Print the likelihood values 
print_likelihood_values(path_to_result)

# Print the fisher correlation matrix
path_fish_cor_mat <- "../../Result/Figures/Model2_nlme/Model2_fish_cor_mat.pdf"
plot_heat_map_fisher(path_fish_cor_mat, path_to_result, print_figure = F)

# Plot the small correlation matrix 
path_fish_cor_mat_small <- "../../Result/Figures/Model2_nlme/Model2_fish_cor_mat_small.pdf"
cor_mat <- plot_heat_map_fisher_small(path_fish_cor_mat_small, path_to_result, print_figure = T)

# Table for the population parameters 
print(create_latex_table_pop_param(path_to_result), include.rownames = F)

# Summary of the EBE:s
print_summary_stat_EBE(str_c(path_to_result, "/"))

# Coefficents of variation 
calc_coeff_of_variation(path_to_result)

# -----------------------------------------------------------------------------------------------------------------
# Model 2 short del 
# -----------------------------------------------------------------------------------------------------------------
# Shrinkage
path_to_result <- "../Monolix/Model2_test/Result"
estimated_sigma <- 0.47996
n_param <- 11
calculate_shrinkage_model(path_to_result, estimated_sigma, n_param)

# Obs vs ipred
path_obs_vs_ipred = "../../Result/Figures/Model2_nlme_short_del/Model2_short_obs_vs_IPRED.pdf"
print_obs_vs_ipred(path_obs_vs_ipred, path_to_result, estimated_sigma, print_figure = T)

# Residuals
path_res_time <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_res_time.pdf"
path_res_response <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_res_response.pdf"
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = T)
print_residuals(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = F, save_time = T)

# Check distribtion 
path_hist_IWRES <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_hist_IWRES.pdf"
path_qq_IWRES <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_qq_IWRES.pdf"
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, print_figure = T)
check_dist_IWRES(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, save_hist = T)

# Check individual parameters 
path_hist_ind <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_hist_ind.pdf"
path_qq_ind <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_qq_ind.pdf"
param_names <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "Suc20", "Glc0")
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            n_col = 6, n_param = 11)
check_dist_individual_param(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names, 
                            save_hist = F, n_col = 6, n_param = 11)

# Check correlation individual parameters 
path_corrplot_ind <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_corrplot_ind.pdf"
path_pair_plot_ind <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_pairs_ind.pdf"
plot_list_corr <- check_correlation_ind_param(path_corrplot_ind, path_pair_plot_ind, path_to_result)
pdf(path_corrplot_ind, width = 9, height = 5); corrplot.mixed(plot_list_corr[[1]]); garbage <- dev.off()
pdf(path_pair_plot_ind, width = 9, height = 5); plot_list_corr[[2]]; garbage <- dev.off()

# Check individual fits 
path_ind_fit <- "../../Result/Figures/Model2_nlme_short_del/Model2_ind_fit_"
plot_individual_fit(path_ind_fit, path_to_result, 74, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 73, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 10, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 26, print_fig = T, save_fig = T)
plot_individual_fit(path_ind_fit, path_to_result, 55, print_fig = T, save_fig = T)

# Print the likelihood values 
print_likelihood_values(path_to_result)

# Print the fisher correlation matrix
path_fish_cor_mat <- "../../Result/Figures/Model2_nlme_short_del/Model2_short_fish_cor_mat.pdf"
plot_heat_map_fisher(path_fish_cor_mat, path_to_result, print_figure = F)

# Plot the small correlation matrix 
path_fish_cor_mat_small <- "../../Result/Figures/Model2_nlme_short_del/Model2_fish_cor_mat_small.pdf"
cor_mat <- plot_heat_map_fisher_small(path_fish_cor_mat_small, path_to_result, print_figure = T)

# Table for the population parameters 
print(create_latex_table_pop_param(path_to_result), include.rownames = F)

# Summary of the EBE:s
print_summary_stat_EBE(str_c(path_to_result, "/"))

# Coefficents of variation 
calc_coeff_of_variation(path_to_result)
