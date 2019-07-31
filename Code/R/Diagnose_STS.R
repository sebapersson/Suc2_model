# This script creates diagnosis plots for the STS result by using the functions in Process_sts_result.R 
source("Process_sts_result.R")


# ------------------------------------------------------------------------------------------------------------------------
# Model 1
# ------------------------------------------------------------------------------------------------------------------------
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_pred_data <- "../Matlab/STS/Model1/Result/Predicted_fits.csv"
path_save <- "../../Result/Figures/Model1_STS/STS_model1_obs_vs_pred.pdf"
plot_obs_vs_pred_STS(path_to_obs_data, path_pred_data, path_save, print_fig = T)

# The residuals 
path_to_residuals <- "../Matlab/STS/Model1/Result/Residuals.csv"
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_save_time <- "../../Result/Figures/Model1_STS/STS_model11_res_time.pdf"
path_save_obs <- "../../Result/Figures/Model1_STS/STS_model11_res_obs.pdf"
plot_residuals_STS(path_to_residuals, path_to_obs_data, path_save_time, path_save_obs, save_time = T, print_fig = T)
plot_residuals_STS(path_to_residuals, path_to_obs_data, path_save_time, path_save_obs, save_time = F, print_fig = T)

# The individual fits 
path_fit_data <- "../Matlab/STS/Model1/Result/Individual_fits.csv"
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_save <- "../../Result/Figures/Model1_STS/STS_model1_ind_fit_"
plot_individual_fit(path_fit_data, path_obs_data, path_save, 74, print_fig = T)
plot_individual_fit(path_fit_data, path_obs_data, path_save, 73, print_fig = T)
plot_individual_fit(path_fit_data, path_obs_data, path_save, 10, print_fig = T)

# The individual parameters
n_parameters <- 12
name_param <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc2", "Glc0", "Cost")
path_param <- "../Matlab/STS/Model1/Result/Parameters_and_cost_func.csv"
path_save_hist <- "../../Result/Figures/Model1_STS/STS_model1_hist_ind_param.pdf"
path_save_qq <- "../../Result/Figures/Model1_STS/STS_model1_qq_ind_param.pdf"
check_dist_ind_param_STS(path_param, name_param, n_parameters, path_save_hist, path_save_qq, save_hist = T)
check_dist_ind_param_STS(path_param, name_param, n_parameters, path_save_hist, path_save_qq, save_hist = F)

# Calculate the MLE estimates for the paremters 
name_param <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc2", "Glc", "Cost")
path_param <- "../Matlab/STS/Model1/Result/Parameters_and_cost_func.csv"
n_parameters <- 12
path_save_file <- "../../Result/Files/STS_model1_cov_mean.csv"
calc_mean_and_cov_param(path_param, name_param, n_parameters, path_save_file)

# ------------------------------------------------------------------------------------------------------------------------
# Model 2
# ------------------------------------------------------------------------------------------------------------------------
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_pred_data <- "../Matlab/STS/Model2/Result/Predicted_fits.csv"
path_save <- "../../Result/Figures/Model2_STS/STS_model2_obs_vs_pred.pdf"
plot_obs_vs_pred_STS(path_to_obs_data, path_pred_data, path_save, print_fig = T)

# The residuals 
path_to_residuals <- "../Matlab/STS/Model2/Result/Residuals.csv"
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_save_time <- "../../Result/Figures/Model2_STS/STS_model2_res_time.pdf"
path_save_obs <- "../../Result/Figures/Model2_STS/STS_model2_res_obs.pdf"
plot_residuals_STS(path_to_residuals, path_to_obs_data, path_save_time, path_save_obs, save_time = T, print_fig = T)
plot_residuals_STS(path_to_residuals, path_to_obs_data, path_save_time, path_save_obs, save_time = F, print_fig = T)

# The individual fits 
path_fit_data <- "../Matlab/STS/Model2/Result/Individual_fits.csv"
path_to_obs_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_save <- "../../Result/Figures/Model2_STS/STS_model2_ind_fit_"
plot_individual_fit(path_fit_data, path_obs_data, path_save, 74, print_fig = T)
plot_individual_fit(path_fit_data, path_obs_data, path_save, 73, print_fig = T)
plot_individual_fit(path_fit_data, path_obs_data, path_save, 10, print_fig = T)

# The individual parameters
n_parameters <- 12
name_param <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc2", "Glc0", "Cost")
path_param <- "../Matlab/STS/Model1/Result/Parameters_and_cost_func.csv"
path_save_hist <- "../../Result/Figures/Model2_STS/STS_model2_hist_ind_param.pdf"
path_save_qq <- "../../Result/Figures/Model2_STS/STS_model2_qq_ind_param.pdf"
check_dist_ind_param_STS(path_param, name_param, n_parameters, path_save_hist, path_save_qq, save_hist = T)
check_dist_ind_param_STS(path_param, name_param, n_parameters, path_save_hist, path_save_qq, save_hist = F)

# Calculate the MLE estimates for the paremters 
name_param <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "Suc2", "Glc", "Cost")
path_param <- "../Matlab/STS/Model2/Result/Parameters_and_cost_func.csv"
n_parameters <- 12
path_save_file <- "../../Result/Files/STS_model2_cov_mean.csv"
calc_mean_and_cov_param(path_param, name_param, n_parameters, path_save_file)
