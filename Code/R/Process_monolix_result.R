library(tidyverse)
library(RColorBrewer)
library(latex2exp)
library(stringr)
library(corrplot)
library(reshape2)
library(xtable)
library(corrr)

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
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# This file contains the functions that take the exported data from monolix and creates graphs from it. 
# Note that this script should be run from the Code/R-folder
# This file contains functions for:


# Function that will calculate, and print, the epsilon and eta shrinkage for a model
# Inputs:
#   path_to_result, the path to the result folder
#   estimated_sigma, the estimated standard deviation for the error model 
#   n_param, the number of parameters in the model 
# Reaturns:
#   void 
calculate_shrinkage_model <- function(path_to_result, estimated_sigma, n_param)
{
  # Calculate the epsilon and eta shrinkage 
  path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/Suc2__obsVsPred.txt")
  obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
    select(-censored, split, color, filter) %>%
    rename("IPRED" = "indivPredMode") %>%
    mutate(IWRES = (Suc2_ - IPRED) / estimated_sigma)
  
  # Calculate the shrinkage, very small which is expected in a data rich setting. 
  eps_shrinkage <- 1 - sd(obs_vs_pred_data$IWRES)
  
  print(sprintf("Epsilon shrinkage = %.3f", eps_shrinkage))
  
  # Calculate the eta shrinkage
  path_to_individual_eta <- str_c(path_to_result, "/IndividualParameters/estimatedRandomEffects.txt")
  individual_eta_raw <- read_csv(path_to_individual_eta, col_types = cols())
  # Only interested in EBE estimates (mode)
  individual_eta_sd <- individual_eta_raw %>% 
    select(ends_with("mode")) %>%
    apply(2, sd)
  # Get the population values for calculations 
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  # Select the omega values  
  pop_val_omega <- pop_val_raw %>%
    filter(grepl("^omega", parameter)) %>%
    select(value)
  # Eta-shrinkage calculation 
  eta_shrinkage <- 1 - individual_eta_sd / (pop_val_omega$value)
  
  # Print the eta-shrinkage
  print("Eta-shrinkage:")
  for(i in 1:n_param){
    k_param <- "k"
    k_param <- str_c(k_param, as.character(i))
    print(sprintf("Shrinkage %s = %.3f", k_param, eta_shrinkage[i]))
  }
}


# Function that will plot the observations vs ipred plot, the figure will then be saved in an apprioate folder
# given by the user, note, the figure will be saved as pdf. The user also gets the alternative to print 
# the figure. 
# Args:
#   path_to_figure, path to where the figure should be saved 
#   path_to_result, the path to the result folder
#   estimated_sigma, the estimated standard deviation for the error model 
#   print_figure, wheter or not the figure should be printed, false by default 
# Returns:
#   void 
print_obs_vs_ipred <- function(path_to_figure, path_to_result, estimated_sigma, print_figure = F)
{
  # Read and process the data 
  path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/Suc2__obsVsPred.txt")
  obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
    select(-censored, split, color, filter) %>%
    rename("IPRED" = "indivPredMode") %>%
    mutate(IWRES = (Suc2_ - IPRED) / estimated_sigma)
  
  # Plotting observation plot 
  cols = c("y = x" = cbPalette[4], "Trend line" = cbPalette[6])
  p1 <- ggplot(obs_vs_pred_data, aes(y = Suc2_, x = IPRED)) + 
    geom_point(alpha = 0.6, size = 0.8) + 
    geom_smooth(aes(color = "Trend line"), method = "lm", size = 1.0, se = F, show.legend = NA) + 
    geom_abline(aes(slope = 1, intercept = 0, color = "y = x"), size = 1.0, show.legend = F) + 
    labs(x = "Individual predictions", y = "Observed values") +
    scale_color_manual(name = "Linetype", values = cols) + 
    my_classic_theme 
  
  # Save the figure 
  p1
  ggsave(path_to_figure, height = 6, width = 9)
  
  # Print the figure if the option provided by the user 
  if(print_figure) print(p1)
  
  return(0)
}


# Function that will plot the IWRES residuals vs time and the response. If given as an option the two 
# graphs will be printed in the same figure. Note that both are written to separate files. 
# Args:
#   path_res_time, the path to the time residuals are saved
#   path_res_response, the path to where the response residuals are saved
#   path_to_result, path to the result folder from monolix 
#   estimated_sigma, the estimated standard deviation for the error model 
#   print_figure, wheter or not the figure should be printed, false by default
#   save_time, if true the residuals vs time will be saved, else the residuals vs response are saved 
# Returns:
#   void 
print_residuals <- function(path_res_time, path_res_response, path_to_result, estimated_sigma, print_figure = F, save_time = F)
{
  # Read and process the data 
  path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/Suc2__obsVsPred.txt")
  obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
    select(-censored, split, color, filter) %>%
    rename("IPRED" = "indivPredMode") %>%
    mutate(IWRES = (Suc2_ - IPRED) / estimated_sigma)
  
  # Plot IWRES vs time 
  p1 <- ggplot(obs_vs_pred_data, aes(time, IWRES)) + 
    geom_point(size = 0.8, alpha = 0.6) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
    geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
    geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
    labs(x = "Time", y = "IWRES") +
    my_classic_theme
  
  # Plot IWRES vs Suc2 expression 
  p2 <- ggplot(obs_vs_pred_data, aes(Suc2_, IWRES)) + 
    geom_point(size = 0.8, alpha = 0.6) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
    geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
    geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
    labs(x = "Suc2-intensity", y = "IWRES") +
    my_classic_theme
  
  # Save the result 
  if(save_time == T){
    # Plot IWRES vs time 
    p1 <- ggplot(obs_vs_pred_data, aes(time, IWRES)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
      labs(x = "Time", y = "IWRES") +
      my_classic_theme
    p1
    ggsave(path_res_time, height = 6, width = 9)
  }else{
    # Plot IWRES vs Suc2 expression 
    p2 <- ggplot(obs_vs_pred_data, aes(Suc2_, IWRES)) + 
      geom_point(size = 0.8, alpha = 0.6) + 
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) + 
      geom_hline(yintercept = -2, linetype = "dashed", size = 0.5) + 
      geom_hline(yintercept = 2, linetype = "dashed", size = 0.5) + 
      geom_smooth(method = "loess", se = F, color = cbPalette[4], size = 1.2) + 
      labs(x = "Suc2-intensity", y = "IWRES") +
      my_classic_theme
    p2 
    ggsave(path_res_response, height = 6, width = 9)
  }
  
  # Print if provided as an option 
  if(print_figure) print(ggpubr::ggarrange(p1, p2, ncol = 2))
  
  return(0)
}


# Function that will create a histogram and qqplot of the IWRES distribution which for a perfect model should 
# be N(0, 1)
# Args:
#   path_hist_IWRES, the path to the histogram is stored 
#   path_qq_IWRES, the path to where the qqplot is stored
#   path_to_result, path to the result folder from monolix 
#   estimated_sigma, the estimated standard deviation for the error model 
#   n_bins, the  number of bins to use for the histogram 
#   print_figure, wheter or not the figure should be printed, false by default 
#   save_hist, if true the histogram will be saved, else the qq-plot is saved 
# Returns:
#   void 
check_dist_IWRES <- function(path_hist_IWRES, path_qq_IWRES, path_to_result, estimated_sigma, n_bins = 100, print_figure = F, save_hist = F)
{
  # Read and process the data 
  path_to_ipred <- str_c(path_to_result, "/ChartsData/ObservationsVsPredictions/Suc2__obsVsPred.txt")
  obs_vs_pred_data <- read_csv(path_to_ipred, col_types = cols()) %>%
    select(-censored, split, color, filter) %>%
    rename("IPRED" = "indivPredMode") %>%
    mutate(IWRES = (Suc2_ - IPRED) / estimated_sigma)
  
  # Simulate the trendline 
  data_seq <- tibble(x_val = seq(from = min(obs_vs_pred_data$IWRES) - 0.1, to = max(obs_vs_pred_data$IWRES) + 0.1, 
                                 by = 0.01)) %>%
    mutate(prob = dnorm(x_val))
  
  # Make a hisotgram of the IWRES values 
  p1 <- ggplot(obs_vs_pred_data, aes(IWRES)) + 
    geom_histogram(aes(y = ..density..), bins = n_bins, position = "identity", color = "black") +
    geom_line(data = data_seq, aes(x_val, prob), color = cbPalette[4], size = 1.2) +
    labs(y = "", x = "IWRES") +
    my_classic_theme
  
  # Create a qqplot 
  p2 <- ggplot(obs_vs_pred_data, aes(sample = IWRES)) + 
    geom_qq() + 
    geom_qq_line() + 
    my_classic_theme
  
  # Save the figures 
  if(save_hist){
    # Make a hisotgram of the IWRES values 
    p1 <- ggplot(obs_vs_pred_data, aes(IWRES)) + 
      geom_histogram(aes(y = ..density..), bins = n_bins, position = "identity", color = "black") +
      geom_line(data = data_seq, aes(x_val, prob), color = cbPalette[4], size = 1.2) +
      labs(y = "", x = "IWRES") +
      my_classic_theme
    p1
    ggsave(path_hist_IWRES, height = 6, width = 9)
  }else{
    # Create a qqplot, this is required, kind of a bug in R 
    p2 <- ggplot(obs_vs_pred_data, aes(sample = IWRES)) + 
      geom_qq() + 
      geom_qq_line() + 
      my_classic_theme
    p2 
    ggsave(path_qq_IWRES, height = 6, width = 9)
  }
  
  # Print if provided as an option 
  if(print_figure) print(ggpubr::ggarrange(p1, p2, ncol = 2))
  
  return(0)
}

# Function that will compute a histogram and qq-plot for the individual parameters. Note that everything 
# is lumped into two big figures. The column and rows for ggarrange have to be provided 
# Args:
#   path_hist_ind, where to save the individual parameters 
#   path_qq_int, where to save the qqplot of the individual parameters. 
#   path_to_result, path to the result folder from monolix 
#   estimated_sigma, the estimated standard deviation for the error model 
#   n_bins, the  number of bins to use for the histogram 
#   n_row, the number of rows for the plot
#   n_col, the number of columns for the plot 
#   Returns:
#   void
check_dist_individual_param <- function(path_hist_ind, path_qq_ind, path_to_result, estimated_sigma, param_names,
                                        n_bins = 50, n_row = 2, n_col = 5, n_param = 10, save_hist = T)
{
  # Read in the sampled parameters 
  path_sample_param <- str_c(path_to_result, "/IndividualParameters/simulatedIndividualParameters.txt")
  sampled_param_raw <- read_csv(path_sample_param, col_types = cols())
  
  # Read in the population values 
  path_pop_val <- str_c(path_to_result, "/populationParameters.txt")
  pop_val_raw <- read_csv(path_pop_val, col_types = cols())
  
  plot_list_hist <- lapply(1:n_param, function(i){
    # Get the sampled values
    sampled_val <- sampled_param_raw %>% 
      select(param_names[i])
    
    # Get population value 
    k_param_str <- str_c(param_names[i], "_pop")
    k_param_i <- which(pop_val_raw$parameter == k_param_str)
    k_val <- pop_val_raw$value[k_param_i]
    
    # Get omega (variance) value
    omega_param_str <- str_c("omega_", param_names[i])
    omega_param_i <- which(pop_val_raw$parameter == omega_param_str)
    omega_val <- pop_val_raw$value[omega_param_i]
    
    # Simulate many normal values
    seq_data <- tibble(x_val = seq(from = min(log(sampled_val)) - 0.1, to = max(log(sampled_val)) + 0.1, by = 0.05)) %>%
      mutate(prob_val = dnorm(x_val, mean = log(k_val), sd = omega_val))
    
    names(sampled_val) <- "to_plot"
    p1 <- ggplot(sampled_val, aes(log(to_plot))) + 
      geom_histogram(aes(y = ..density..), bins = n_bins, position = "identity", color = "black") +
      geom_line(data = seq_data, aes((x_val), prob_val), color = cbPalette[4], size = 1.2) + 
      labs(y = "", x = param_names[i]) + 
      my_classic_theme
    
    return(p1)})
  
  # Compute the qq-plots 
  plot_list_qq <- lapply(1:n_param, function(i){
    # Get the sampled values
    sampled_val <- sampled_param_raw %>% 
      select(param_names[i])
    
    names(sampled_val) <- "to_plot"
    p1 <- ggplot(sampled_val, aes(sample = log(to_plot))) + 
      geom_qq_line() + 
      geom_qq() + 
      labs(x = param_names[i]) +
      my_classic_theme
    
    return(p1)})
  
  # Plot the histogram, can only plot one...
  if(save_hist){
    p1 <- ggpubr::ggarrange(plotlist = plot_list_hist, ncol = n_col, nrow = 2)
    p1
    ggsave(path_hist_ind, height = 6, width = 9)
  }else{
    # Plot the qqplots
    p2 <- ggpubr::ggarrange(plotlist = plot_list_qq, ncol = n_col, nrow = 2)
    p2 
    ggsave(path_qq_ind, height = 6, width = 9)
  }
  
  return(0)
}

# Function that will compuate the correlation between the individual parameters, and then plot it using pairs 
# and corrplot. 
# Args:
#   path_corrplot_ind, where to save the correlation plots 
#   path_pair_plot_ind, where to save the pairs plot
#   path_to_result, path to the result folder from monolix 
check_correlation_ind_param <- function(path_corrplot_ind, path_pair_plot_ind, path_to_result)
{
  
  # Read the conditional parameters 
  path_sample_param <- str_c(path_to_result, "/IndividualParameters/simulatedIndividualParameters.txt")
  sampled_param_raw <- read_csv(path_sample_param, col_types = cols())
  sampled_param <- sampled_param_raw %>%
    select(-rep, -id)
  
  # Calculate the covariance matrix 
  cor_mat <- cov2cor(cov(sampled_param))
  p1 <- corrplot.mixed(cor_mat)
  
  # Scale the sampled paramaters
  sampled_param_scaled <- as_tibble(as.data.frame(scale(sampled_param)))
  
  # Plot the data
  lowerFn <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point() +
      geom_smooth(method = method, color = cbPalette[4], se = F)
    p
  }
  
  p2 <- GGally::ggpairs(sampled_param_scaled, 
                        lower = list(continuous = GGally::wrap(lowerFn))) + 
    my_classic_theme

  # Save the result
  return(list(p1, p2, cor_mat))
}  

# Function that will plot the predictions (based on EBE) vs observed data for an individual cell. The picture
# can be saved if the user provides that option. 
# Args:
#   path_fig, where to store the result, note which cell is appended to path_fig
#   path_to_result, path to the result folder from monolix 
#   which_cell, which cell to plot 
#   save_fig, if the figure should be saved or not, deafult false
#   print_fig, if the figure should be printed or not, deafult false
# Returns:
#   void
plot_individual_fit <- function(path_fig, path_to_result, which_cell, save_fig = F, print_fig = F)
{
  # Individual fits 
  path_fits <- str_c(path_to_result, "/ChartsData/IndividualFits/Suc2__fits.txt")
  fit_data <- read_csv(path_fits, col_types = cols())
  
  # Observations 
  path_observations <- str_c(path_to_result, "/ChartsData/IndividualFits/Suc2__observations.txt")
  obs_data <- read_csv(path_observations, col_types = cols())
  
  # Name for saving the figure 
  name_fig <- str_c(path_fig, as.character(which_cell), ".pdf")
  
  # Filter out specific cell to plot 
  if(names(obs_data)[1] == "ID"){
    data_fit_cell <- fit_data %>% filter(ID == which_cell)
    data_obs_cell <- obs_data %>% filter(ID == which_cell)
  }else{
    data_fit_cell <- fit_data %>% filter(Cell == which_cell)
    data_obs_cell <- obs_data %>% filter(Cell == which_cell)
  }
  
  # Plot the data 
  p1 <- ggplot(data_fit_cell, aes(time, indivPredMode)) + 
    geom_line(size = 1.1, color = cbPalette[3]) + 
    geom_point(data = data_obs_cell, aes(time, Suc2_), color = cbPalette[2], size = 1.6) + 
    labs(x = TeX("Scaled time"), y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) + 
    my_classic_theme
  
  # Save fig if required 
  if(save_fig){
    p1 
    ggsave(name_fig, height = 6, width = 9)
  }
  
  # If the figure should be printed
  if(print_fig) print(p1)
  
  return(0)
}


# Function that will print the likelihood values for a model, more specifically it will print the likelhood,  
# AIC and BIC. 
# Args:
#   path_to_result, path to the result folder 
# Returns:
#   void
print_likelihood_values <- function(path_to_result)
{
  # Process the likelihood result 
  path_likelihood <- str_c(path_to_result, "/LogLikelihood/logLikelihood.txt")
  likelihood_data <- read_csv(path_likelihood, col_types = cols())
  
  # Print the data 
  print(sprintf("Log-likelhood = %.2f", likelihood_data[1, 2]))
  print(sprintf("AIC = %.2f", likelihood_data[2, 2]))
  print(sprintf("BIC = %.2f", likelihood_data[4, 2]))
  return(0)
}

# Function that will plot a heatmap of the fisher correlation matrix for a model. 
# Args:
#   path_fish_cor, where to save the correlation matrix 
#   path_to_result, the path to where the monolix result is stored. 
#   print_figure, if the figure should be printed or not, false by default 
# Returns:
#   void 
plot_heat_map_fisher <- function(path_fish_cor, path_to_result, print_figure = F)
{

  # Read and process data into correlation matrix 
  path_corr_mat <- str_c(path_to_result, "/FisherInformation/correlationEstimatesSA.txt")
  cor_mat <- read_csv(path_corr_mat, col_types = cols(), col_names = F)
  names(cor_mat)[2:dim(cor_mat)[2]] <- cor_mat$X1
  cor_mat <- cor_mat %>% select(-X1)
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  lower_tri <- get_lower_tri(as.matrix(cor_mat))
  rownames(lower_tri) <- colnames(lower_tri)
  meleted_cormat <- melt(lower_tri, na.rm = T)

  # Print the heat-map
  p1 <- ggplot(data = meleted_cormat, aes(Var1, Var2, fill = value)) + 
    geom_tile(color = "white") + 
    scale_fill_gradient2(low = "#7b3294", high = "#008837", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
    coord_fixed() + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.6),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.0,
                                 title.position = "top", title.hjust = 0.3))
  
  # Save the figure 
  p1
  ggsave(path_fish_cor, height = 6, width = 9)
  
  # If the figure should be printed
  if(print_figure) print(p1)
  
  return(0)
}

# Function that will plot a heatmap of the fisher correlation matrix for a model, but it will only plot the 
# correlation parameters. 
# Args:
#   path_fish_cor, where to save the correlation matrix 
#   path_to_result, the path to where the monolix result is stored. 
#   print_figure, if the figure should be printed or not, false by default 
# Returns:
#   void 
plot_heat_map_fisher_small <- function(path_fish_cor, path_to_result, print_figure = F)
{
  # Read and process data into correlation matrix 
  path_corr_mat <- str_c(path_to_result, "/FisherInformation/correlationEstimatesSA.txt")
  cor_mat <- read_csv(path_corr_mat, col_types = cols(), col_names = F) %>%
    filter(str_detect(X1, "^k"))
  names(cor_mat)[2:(dim(cor_mat)[1] + 1)] <- cor_mat$X1
  print("Got here")
  cor_mat <- cor_mat %>% select(-X1) %>% select(starts_with("k"))
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  lower_tri <- get_lower_tri(as.matrix(cor_mat))
  rownames(lower_tri) <- colnames(lower_tri)
  meleted_cormat <- melt(lower_tri, na.rm = T)
  
  # Print the heat-map
  p1 <- ggplot(data = meleted_cormat, aes(Var1, Var2, fill = value)) + 
    geom_tile(color = "white") + 
    scale_fill_gradient2(low = "#7b3294", high = "#008837", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+
    coord_fixed() + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.6),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 8, barheight = 1.0,
                                 title.position = "top", title.hjust = 0.3))
  
  # Save the figure 
  p1
  ggsave(path_fish_cor, height = 6, width = 9)
  
  # If the figure should be printed
  if(print_figure) print(p1)
  
  return(cor_mat)
}

# Function that will create code for a LaTeX table out of the population parameters. 
# Args:
#   path_to_result, the path to the result folder 
# Returns:
#   void 
create_latex_table_pop_param <- function(path_to_result)
{
  # Read in the population parameters and make latex-table 
  path_pop_param <- str_c(path_to_result, "/populationParameters.txt")
  pop_param <- read_csv(path_pop_param, col_types = cols()) %>%
    mutate(value = as.character(round(value, 2))) %>%
    mutate(rse_sa = as.character(round(rse_sa, 1))) %>%
    select(-se_sa) %>%
    rename("Parameter" = "parameter") %>%
    mutate(value = str_c(value, " (", rse_sa, ")")) %>%
    select(-rse_sa) %>%
    rename("Value" = "value")
  
  # Return the LaTex code (will only have to be altered a little in LaTex)
  latex_table <- xtable(pop_param)
  return(latex_table)
}


# Function that prints the coefficent of variation for a certain model, note that each parameter is assumed to be 
# log normaly distributed. 
# Args:
#   path_to_result, the path to the result folder 
# Returns:
#   void 
calc_coeff_of_variation <- function(path_to_result)
{
  # Calculate the coefficent of variation and output the result 
  path_pop_param <- str_c(path_to_result, "/populationParameters.txt")
  pop_data <- read_csv(path_pop_param, col_types = cols())
  omega_data <- pop_data %>%
    filter(grepl("omega", parameter)) %>%
    select(parameter, value) %>%
    mutate(Coeff_var = sqrt(exp(value^2) - 1))
  
  # Get the names of the parmeters 
  parameter_names <- sub("^[^_]*", "", omega_data$parameter) %>%
    substring(2)
  
  # Print the result 
  for(i in 1:length(parameter_names)){
    print(sprintf("%s: %.3f", parameter_names[i], omega_data$Coeff_var[i]))
  }
  
  return(0)
}

# Function that will calculate the correlation between the EBE:s and the area of the cells (the size with other words). 
# Args:
#   path_to_result, the path to the monolix result folder
# Returns:
#   void 
calc_corr_area <- function(path_to_result)
{

  # Read the EBE:s 
  path_param <- str_c(path_to_result, "IndividualParameters/estimatedIndividualParameters.txt")
  param_data_raw <- read_csv(path_param, col_types = cols(id = col_factor())) %>%
    select(id, ends_with("mode"))
  
  # Read the entire data-set 
  data_tot <- read_csv("../../Intermediate/Data_whole_tidy_filt_area.csv", col_types = cols()) %>%
    filter(t == 0) %>%
    select(Area)
  
  # Add area to the EBE:s
  param_area_data <- param_data_raw %>%
    bind_cols(data_tot) %>% 
    select(-id)
  
  # Correlation data 
  corr_result <- param_area_data %>% 
    correlate(method = "pearson", quiet = T) %>%
    focus(Area) %>%
    rename("Parameter" = "rowname")
  
  # Print the final result in LaTeX table format 
  print(xtable(t(corr_result)))
  
  return(0)
}


# Function that will print a LaTex table of summarise statistic (mean, meadian, min, max, 10 % and 90 % quantile)
# of the EBE:s for a model.
# Args:
#   path_to_result, the path to the monolix result folder
# Returns:
#   void 
print_summary_stat_EBE <- function(path_to_result)
{
  # Read the EBE:s 
  path_param <- str_c(path_to_result, "IndividualParameters/estimatedIndividualParameters.txt")
  param_data_raw <- read_csv(path_param, col_types = cols(id = col_factor())) %>%
    select(id, ends_with("mode"))
  
  # Change to long format for calculation 
  param_data_tidy <- param_data_raw %>%
    gather(ends_with("mode"), value = "Param_val", key = "Parameter") %>%
    select(-id)
  
  # Calculate the summarise statistics 
  param_data_sum <- param_data_tidy %>%
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

