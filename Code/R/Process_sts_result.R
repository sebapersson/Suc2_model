library(tidyverse)
library(RColorBrewer)
library(latex2exp)
library(stringr)
library(corrplot)
#library(MASS)
#library(ggpubr)

RNGkind("L'Ecuyer-CMRG")

# General plotting parameters 
my_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                               plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13))
my_min_theme <- theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13))
my_classic_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                            plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13)) + theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                                                                  colour = "grey"))

# For certain kinds of plots where one parameter increses 
my_palette <- brewer.pal(n = 9, name = "GnBu")[-c(1, 2, 3)]

# colour-blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# This file contains functions for processing the STS result. The function will store the created functions, files etc 
# in a location provided by the user. 

# Function that will plot the observed vs fitted values for all the individual fits. The plot will have an identity line
# and a trend line. Note, the datais assumed to have a certain standard format. 
# Args:
#   path_obs_data, path to the observed data 
#   path_pred_data, path to the predicted fits for the model 
#   path_save, where to store the figure 
#   print_fig, if true the figure is ploted, by default false 
# Returns:
#   void
plot_obs_vs_pred_STS <- function(path_to_obs_data, path_pred_data, path_save, print_fig = F)
{
  # Read and fix data-types for the data 
  obs_data <- read_csv(path_to_obs_data, col_types = cols(ID = col_factor())) %>%
    select(t, Suc2, ID) 
  pred_data <- read_csv(path_pred_data, col_names = F, col_types = cols()) %>%
    rename("t" = "X1", "Suc2_pred" = "X2", "ID" = "X3") %>%
    mutate(ID = as.factor(ID))
  
  # Create tible to plot with observed and predicted
  data_to_plot <- obs_data %>%
    select(Suc2) %>%
    bind_cols(pred_data) %>%
    select(t, Suc2, Suc2_pred, ID)
  
  # Draw the graph   
  cols = c("y = x" = cbPalette[4], "Trend line" = cbPalette[6])
  p1 <- ggplot(data_to_plot, aes(y = Suc2, x = Suc2_pred)) + 
    geom_point(alpha = 0.6, size = 0.8) + 
    geom_smooth(aes(color = "Trend line"), method = "lm", size = 1.0, se = F, show.legend = NA) + 
    geom_abline(aes(slope = 1, intercept = 0, color = "y = x"), size = 1.0, show.legend = F) + 
    labs(x = "Individual predictions", y = "Observed values") +
    scale_color_manual(name = "Linetype", values = cols) + 
    my_classic_theme 
  
  # Save figures
  p1 
  ggsave(path_save, height = 6, width = 9)
  
  # Print if provded as option 
  if(print_fig) print(p1)
  
  return(0) 
}


# Function that will plot all residuals vs time and vs the observed values. Due to technical reasons in R only 
# one of the plots can be printed and saved by call. 
# Args:
#   path_to_residuals, the path to the residuals 
#   path_obs_data, path to the observed data 
#   path_save_time, the path where to save the time residuals
#   path_save_obs, the path where to save the observed residuals 
#   save_time, if the time residuals should be saved, by default T (this reuslt in obs not being saved)
#   print_fig, if the figure should be printed or not, by default false 
# Returns:
#   void 
plot_residuals_STS <- function(path_to_residuals, path_obs_data, path_save_time, path_save_obs, 
                               save_time = T, print_fig = F)
{
  # Read and process the data 
  residuals_data <- read_csv(path_to_residuals, col_names = F, col_types = cols()) %>%
    gather(starts_with("X"), value = Residual, key = ID_res) %>%
    select(-ID_res)
  obs_data <- read_csv(path_to_obs_data, col_types = cols(ID = col_factor())) %>%
    select(t, Suc2, ID) 
  # Store the plotting data in a tibble 
  residual_data_plot <- obs_data %>%
    bind_cols(residuals_data) %>%
    mutate(std_res = Residual / sd(Residual))
  
  # Save the choosen figure 
  if(save_time){
    p1 <- ggplot(residual_data_plot, aes(t, std_res)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
      labs(x = "Time", y = "Stand. residuals") +
      my_classic_theme
    p1 
    ggsave(path_save_time, height = 6, width = 9)
  }else{
    p1 <- ggplot(residual_data_plot, aes(Suc2, std_res)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
      labs(x = "Time", y = "Stand. residuals") +
      my_classic_theme
    p1
    ggsave(path_save_obs, height = 6, width = 9)
  }
  
  # Print the figure if provided as an option 
  if(print_fig) print(p1)
  
  return(0)
}


# Function that will plot the inidividual fit for a cell and then save the figure. Note that all data files are 
# assumed to follow a standard format. 
# Args:
#   path_fit_data, the path to the computed fits 
#   path_obs_data, the path to the observed data 
#   path_save, the path where to save the file (which cell is appended to the name)
#   which_cell, the index of the cell to plot, has to be below the maximym number of cells 
#   print_fig, if the figure should be printed or not, false by default
# Returns:
#   void
plot_individual_fit <- function(path_fit_data, path_obs_data, path_save, which_cell, print_fig = F)
{
  # Read and process the data, assumed to be in a standard format 
  fit_data <- read_csv(path_fit_data, col_names = F, col_types = cols()) %>%
    rename("t" = "X1", "Suc2_fit" = "X2", "ID" = "X3") %>%
    mutate(ID = as.factor(ID)) %>%
    filter(ID == which_cell)
  obs_data <- read_csv(path_to_obs_data, col_types = cols(ID = col_factor())) %>%
    select(t, Suc2, ID) %>%
    filter(ID == which_cell)
  title_text <- str_c("Individual fit for cell ", as.character(which_cell))
  
  p1 <- ggplot(obs_data, aes(t, Suc2)) + 
    geom_point(color = cbPalette[2], size = 1.6) + 
    geom_line(fit_data, mapping = aes(t, Suc2_fit), color = cbPalette[3], size = 1.1) + 
    labs(x = "Scale time", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), title = title_text) + 
    my_classic_theme
  
  # Save the figure 
  path_save_fig <- str_c(path_save, as.character(which_cell), ".pdf")
  p1
  ggsave(path_save_fig, height = 6, width = 9)
  
  # Print the figure if provided as option 
  if(print_fig) print(p1)
  
  return(0)
}


# Function that will plot a histogram (with a fitted log-normal) and qq-plots for each individual parameter. 
# Due to the workings of R only one of the plots can be saved in each run. 
# Args:
#   path_param, path to the paremeters file 
#   name_param, name for the parameters, this mainly differes between model 6 and 7, last column should be Cost!
#   n_parameters, the  number of parameters in the model 
#   path_save_hist, the path to where to save the histograms
#   path_save_qq, the path to where to save the qq-plot 
#   save_hist, if to save the histogram or not, true by default
# Returns:
#   void
check_dist_ind_param_STS <- function(path_param, name_param, n_parameters, path_save_hist, path_save_qq, save_hist = T)
{
  # Read and proces the data
  param_data <- read_csv(path_param, col_names = F, col_types = cols())
  names(param_data) <- name_param
  param_data <- param_data %>% select(-Cost)
  
  plot_list_hist <- lapply(1:n_parameters, function(i){
    # For being able to use ggplot 
    data_plot <- param_data[, i]
    names(data_plot) <- "data"
    
    # Fit a log-normal distribution 
    fit_est <- MASS::fitdistr(data_plot$data, "lognormal")$estimate
    dist_data_log <- tibble(x = seq(from = min(data_plot$data) - 0.1, to = max(data_plot$data) + 0.1, by = 0.05)) %>%
      mutate(prob = dlnorm(x, meanlog = fit_est[1], sdlog = fit_est[2]))
    
    p1 <- ggplot(data_plot, aes(data)) + 
      geom_histogram(aes(y = ..density..), bins = 30, color = "black") + 
      geom_line(dist_data_log, mapping = aes(x, prob), color = cbPalette[4], size = 1.2) + 
      labs(y = "Density", x = str_c("Parameter: ", name_param[i])) + 
      my_classic_theme
    return(p1)})
  
  plot_list_qq <- lapply(1:n_parameters, function(i){
    # For being able to use ggplot, log is taken in order to create a qq-plot  
    data_plot <- log(param_data[, i])
    names(data_plot) <- "data"
    
    p1 <- ggplot(data_plot, aes(sample = data)) + 
      geom_qq() + 
      geom_qq_line()+ 
      labs(title = str_c("Parameter: ", name_param[i]), "qq plot") + 
      my_classic_theme
    return(p1)})
  
  if(save_hist){
    p1 <- ggpubr::ggarrange(plotlist = plot_list_hist, nrow = 2, ncol = 6)
    p1
    ggsave(path_save_hist, height = 5, width = 12)
    return(plot_list_hist)
  }else{
    p1 <- ggpubr::ggarrange(plotlist = plot_list_qq, nrow = 2, ncol = 6)
    p1
    ggsave(path_save_qq, height = 5, width = 12)
    return(plot_list_qq)
  }
  return(0)
}


# Function that will compuate the mean vector and sample covariance matrix for the logarithm of the random effects. The resulting 
# parameters will then be written to file. 
# Args:
#   path_param, the path to the calculates parameters 
#   name_param, the name of the parameters
#   n_parameters, the number of parameters in the model
#   path_save_file, where to store the resulting mean vector and covariance matrix 
# Returns:
#   void 
calc_mean_and_cov_param <- function(path_param, name_param, n_parameters, path_save_file)
{
  # Read and process data 
  param_data <- read_csv(path_param, col_names = F, col_types = cols())
  names(param_data) <- name_param
  param_data <- param_data %>% select(-Cost)
  
  # Turn into matrix for compuations and log each value to work with normaly distributed data 
  param_data_mat <- log(as.matrix(param_data))
  
  # Compute the number of observations 
  n_obs <- dim(param_data_mat)[1]
  
  # Computing the mean vector
  mean_vec <- colMeans(param_data_mat)
  
  # Compuate the covariance matrix sigma, will use the sample covariance matrix to get unbiased estimate
  sigma_hat <- matrix(0, ncol = n_parameters, nrow = n_parameters)
  for(i in 1:n_obs){ 
    sigma_hat <- sigma_hat + (param_data_mat[i, ] - mean_vec) %*% t(param_data_mat[i, ] - mean_vec)
  }
  sigma_hat <- sigma_hat / (n_obs - 1)
  
  # Write the data to disk 
  data_save <- as_tibble(sigma_hat) %>%
    bind_rows(mean_vec)
  write_csv(data_save, path_save_file)
  
  return(0)
}


# Function that will print summary statistics (mean, meadian, min, max, 10 % and 90 % quantile) for the STS
# estimated parameters. 
# Args:
#   path_param, path to the parameters folder
#   names_param, name for the parameters in correct order 
# Returns:
#   void 
calc_summary_stats_STS_param <- function(path_param, name_param)
{
  # Read and process the data into long format 
  param_data <- read_csv(path_param, col_names = F, col_types = cols())
  names(param_data) <- names_param
  param_data <- param_data %>%
    select(-Cost) %>%
    gather(value = "Param_val", key = "Parameter")
  
  # Calculate the summarise statistics 
  param_data_sum <- param_data %>%
    group_by(Parameter) %>%
    summarise(Mean = mean(Param_val), 
              Median = median(Param_val), 
              Min = min(Param_val), 
              Max = max(Param_val), 
              Quant10 = quantile(Param_val, 0.1), 
              Quant90 = quantile(Param_val, 0.9))
  
  # Create LaTex code of the table
  print(xtable(param_data_sum))
  
  return(0)
}