library(tidyverse)
library(latex2exp)
library(RColorBrewer)
library(stringr)

RNGkind("L'Ecuyer-CMRG")
my_classic_theme <- theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                            plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13)) + theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                                                                  colour = "grey"))
my_palette <- brewer.pal(n = 9, name = "Blues")[-c(1, 2, 3)]

# This file will plot the result from running the simulation of cells. It will produce a plot for each state in the 
# model. The choosen states can then be savade

# Function that will plot the simulated data for a four state model. The simulated data will be plotted in form 
# of two quantiles provded by the user. For the outsignal Suc2 the function will plot the observed data together with 
# the simulated data, and it will plot the observed utermost quantiles and median. 
# Args:
#   path_to_data, path to the simulated data 
#   state_names, for for the different states in the model
# Returns:
#   a list with the different four plots, the plots in the list can then be exported 
plot_simulated_extrapolation_data <- function(path_to_data, state_names)
{
  # Quantiles to plot 
  quant1 <- c(0.20, 0.80)
  quant2 <- c(0.05, 0.95)
  
  # Read and process the data simulated data
  simulated_data <- read_csv(path_to_data, col_names = F, col_types = cols()) 
  names(simulated_data) <- state_names
  simulated_data <- simulated_data %>% mutate(Cell = as.factor(Cell))
  
  # Summarise the data 
  simulated_data_sum <- simulated_data %>%
    group_by(t) %>%
    summarise(Glc_low1 = quantile(Glc, quant1[1], na.rm = T), 
              Glc_low2 = quantile(Glc, quant2[1], na.rm = T), 
              Glc_high1 = quantile(Glc, quant1[2], na.rm = T), 
              Glc_high2 = quantile(Glc, quant2[2], na.rm = T),
              Glc_median = median(Glc, na.rm = T), 
              Mig1_low1 = quantile(Mig1, quant1[1], na.rm = T), 
              Mig1_low2 = quantile(Mig1, quant2[1], na.rm = T), 
              Mig1_high1 = quantile(Mig1, quant1[2], na.rm = T), 
              Mig1_high2 = quantile(Mig1, quant2[2], na.rm = T), 
              Mig1_median = median(Mig1, na.rm = T), 
              Suc2_low1 = quantile(Suc2, quant1[1], na.rm = T), 
              Suc2_low2 = quantile(Suc2, quant2[1], na.rm = T), 
              Suc2_high1 = quantile(Suc2, quant1[2], na.rm = T), 
              Suc2_high2 = quantile(Suc2, quant2[2], na.rm = T), 
              Suc2_median = median(Suc2, na.rm = T), 
              X_low1 = quantile(X, quant1[1], na.rm = T), 
              X_low2 = quantile(X, quant2[1], na.rm = T), 
              X_high1 = quantile(X, quant1[2], na.rm = T), 
              X_high2 = quantile(X, quant2[2], na.rm = T), 
              X_median = median(X, na.rm = T))
  
  # Plot the glucose (or first state)
  p1 <- ggplot(simulated_data_sum, aes(t, Glc_median)) + 
    geom_line(color = my_palette[6], size = .8) +
    geom_ribbon(aes(ymin = Glc_low1, ymax = Glc_high1), alpha = 0.6, fill = my_palette[4], color = NA) + 
    geom_ribbon(aes(ymin = Glc_high1, ymax = Glc_high2), alpha = 0.6, fill = my_palette[2], color = NA) + 
    geom_ribbon(aes(ymin = Glc_low2, ymax = Glc_low1), alpha = 0.6, fill = my_palette[2], color = NA) + 
    labs(y = str_c(state_names[1]), 
         title = str_c(state_names[1], " vs time"),  x = "Scaled time") + 
    my_classic_theme
  
  # Ploting state 2
  p2 <- ggplot(simulated_data_sum, aes(t, Mig1_median)) + 
    geom_line(color = my_palette[6], size = .8) +
    geom_ribbon(aes(ymin = Mig1_low1, ymax = Mig1_high1), alpha = 0.6, fill = my_palette[4], color = NA) + 
    geom_ribbon(aes(ymin = Mig1_high1, ymax = Mig1_high2), alpha = 0.6, fill = my_palette[2], color = NA) + 
    geom_ribbon(aes(ymin = Mig1_low2, ymax = Mig1_low1), alpha = 0.6, fill = my_palette[2], color = NA) + 
    labs(y = str_c(state_names[3]), 
         title = str_c(state_names[3], " vs time"),  x = "Scaled time") + 
    my_classic_theme
  
  # Plotting state 3
  p3 <- ggplot(simulated_data_sum, aes(t, Suc2_median)) + 
    geom_line(color = my_palette[6], size = .8) +
    geom_ribbon(aes(ymin = Mig1_low1, ymax = Suc2_high1), alpha = 0.6, fill = my_palette[4], color = NA) + 
    geom_ribbon(aes(ymin = Mig1_high1, ymax = Suc2_high2), alpha = 0.6, fill = my_palette[2], color = NA) + 
    geom_ribbon(aes(ymin = Mig1_low2, ymax = Suc2_low1), alpha = 0.6, fill = my_palette[2], color = NA) + 
    labs(y = str_c(state_names[2]), 
         title = str_c(state_names[2], " vs time"),  x = "Scaled time") + 
    my_classic_theme
  
  p4 <- ggplot(simulated_data_sum, aes(t, X_median)) + 
    geom_line(color = my_palette[6], size = .8) +
    geom_ribbon(aes(ymin = X_low1, ymax = X_high1), alpha = 0.6, fill = my_palette[4], color = NA) + 
    geom_ribbon(aes(ymin = X_high1, ymax = X_high2), alpha = 0.6, fill = my_palette[2], color = NA) + 
    geom_ribbon(aes(ymin = X_low2, ymax = X_low1), alpha = 0.6, fill = my_palette[2], color = NA) + 
    labs(y = str_c(state_names[4]), 
         title = str_c(state_names[4], " vs time"),  x = "Scaled time") + 
    my_classic_theme
  
  plot_list <- list(p1, p2, p3, p4) 
  
  return(plot_list)
} 


# Model2 extend time 
state_names <- c("Glc", "Mig1", "Suc2", "X", "t", "Cell")
path_to_data <- "../../Intermediate/Simulation_model2_extrapolate.csv"
plot_list <- plot_simulated_extrapolation_data(path_to_data, state_names)
plot_list[[1]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Glc_time.pdf")
plot_list[[2]] + labs(title = "SNF1 pathway vs time", y = "SNF1 pathway")
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_SNF1_time.pdf")
plot_list[[3]] + labs(y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"))
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Suc2_time.pdf")
plot_list[[4]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_X_time.pdf")

# Model2 rate in 1 / 40 -> 1/ 20
state_names <- c("Glc", "Mig1", "Suc2", "X", "t", "Cell")
path_to_data <- "../../Intermediate/Simulation_model2_extra_glc_1_20.csv"
plot_list <- plot_simulated_extrapolation_data(path_to_data, state_names)
plot_list[[1]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Glc_1_20.pdf")
plot_list[[2]] + labs(title = "SNF1 pathway vs time", y = "SNF1 pathway")
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_SNF1_1_20.pdf")
plot_list[[3]] + labs(y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"))
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Suc2_1_20.pdf")
plot_list[[4]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_X_1_20.pdf")

# Model2 rate in 1 / 40 -> 1/ 20
state_names <- c("Glc", "Mig1", "Suc2", "X", "t", "Cell")
path_to_data <- "../../Intermediate/Simulation_model2_extra_glc_1_2.csv"
plot_list <- plot_simulated_extrapolation_data(path_to_data, state_names)
plot_list[[1]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Glc_1_2.pdf")
plot_list[[2]] + labs(title = "SNF1 pathway vs time", y = "SNF1 pathway")
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_SNF1_1_2.pdf")
plot_list[[3]] + labs(y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"))
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_Suc2_1_2.pdf")
plot_list[[4]]
ggsave("../../Result/Figures/Model2_extrapolated/Simulation_nlme2_X_1_2.pdf")