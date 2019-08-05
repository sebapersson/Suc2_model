library(tidyverse)
library(latex2exp)
library(RColorBrewer)

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

my_palette <- brewer.pal(n = 9, name = "GnBu")[-c(1, 2, 3)]

# Function that will create plots for the VPC test for a certain model. It will provide three plots that 
# will check the (0.05, 0.95), (0.2, 0.8) and (0.3, 0.7) quantiles. 
# Args:
#   path_to_data, the path to the observed data 
#   path_to_vpc_data, the path to the simulated vpc data 
# Returns
#   A list of plots with the intervalls for the simulated quantiles plotted against the observed
create_vpc_plots <- function(path_to_data, path_to_vpc_data)
{
  # Read the raw data and process the quantiles 
  data <- read_csv(path_to_data, col_types = cols()) %>%
    select(t, Suc2)
  data_sum <- data %>%
    group_by(t) %>%
    summarise(median = median(Suc2), 
              quant01 = quantile(Suc2, 0.025), 
              quant02 = quantile(Suc2, 0.2), 
              quant03 = quantile(Suc2, 0.3), 
              quant07 = quantile(Suc2, 0.7), 
              quant08 = quantile(Suc2, 0.8),
              quant09 = quantile(Suc2, 0.975))
  
  # Read the vpc-data 
  vpc_data <- read_csv(path_to_vpc_data, col_types = cols(), col_names = F) %>%
    rename("t" = "X1", "median" = "X2", "quant01" = "X3", "quant02" = "X4", "quant03" = "X5", "quant07" = "X6", 
           "quant08" = "X7", "quant09" = "X8")
  vpc_data_sum <- vpc_data %>%
    group_by(t) %>%
    summarise(upp_quant_median = quantile(median, 0.95, na.rm = T), 
              low_quant_median = quantile(median, 0.05, na.rm = T), 
              low_quant01 = quantile(quant01, 0.05, na.rm = T), 
              upp_quant01 = quantile(quant01, 0.95, na.rm = T), 
              low_quant02 = quantile(quant02, 0.05, na.rm = T), 
              upp_quant02 = quantile(quant02, 0.95, na.rm = T),
              low_quant03 = quantile(quant03, 0.05, na.rm = T), 
              upp_quant03 = quantile(quant03, 0.95, na.rm = T), 
              low_quant07 = quantile(quant07, 0.05, na.rm = T), 
              upp_quant07 = quantile(quant07, 0.95, na.rm = T),
              low_quant08 = quantile(quant08, 0.05, na.rm = T), 
              upp_quant08 = quantile(quant08, 0.95, na.rm = T),
              low_quant09 = quantile(quant09, 0.05, na.rm = T), 
              upp_quant09 = quantile(quant09, 0.95, na.rm = T))
  
  # Create the graphs with confidence bands 
  p1 <- ggplot() + 
    geom_line(data = data_sum, aes(t, median), color = my_palette[6]) + 
    geom_line(data = data_sum, aes(t, quant03), color = my_palette[3]) + 
    geom_line(data = data_sum, aes(t, quant07), color = my_palette[3]) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant03, ymax = upp_quant03), 
                fill = my_palette[3], alpha = 0.3, color = NA) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant07, ymax = upp_quant07), 
                fill = my_palette[3], alpha = 0.3, color = NA) + 
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant_median, ymax = upp_quant_median), 
                fill = my_palette[6], alpha = 0.3, color = NA) + 
    labs(title = "0.3 and 0.7 quantiles", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), 
         x = "Scaled time") + 
    my_classic_theme
  
  p2 <- ggplot() + 
    geom_line(data = data_sum, aes(t, median), color = my_palette[6]) + 
    geom_line(data = data_sum, aes(t, quant02), color = my_palette[3]) + 
    geom_line(data = data_sum, aes(t, quant08), color = my_palette[3]) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant02, ymax = upp_quant02), 
                fill = my_palette[3], alpha = 0.3, color = NA) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant08, ymax = upp_quant08), 
                fill = my_palette[3], alpha = 0.3, color = NA) + 
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant_median, ymax = upp_quant_median), 
                fill = my_palette[6], alpha = 0.3, color = NA) + 
    labs(title = "0.2 and 0.8 quantiles", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), 
         x = "Scaled time") + 
    my_classic_theme
  
  p3 <- ggplot() + 
    geom_line(data = data_sum, aes(t, median), color = my_palette[6]) + 
    geom_line(data = data_sum, aes(t, quant01), color = my_palette[3]) + 
    geom_line(data = data_sum, aes(t, quant09), color = my_palette[3]) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant01, ymax = upp_quant01), 
                fill = my_palette[3], alpha = 0.3, color = NA) +
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant09, ymax = upp_quant09), 
                fill = my_palette[3], alpha = 0.3, color = NA) + 
    geom_ribbon(data = vpc_data_sum, mapping = aes(x = t, ymin = low_quant_median, ymax = upp_quant_median), 
                fill = my_palette[6], alpha = 0.3, color = NA) + 
    labs(title = "0.025 and 0.975 quantiles", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), 
         x = "Scaled time") + 
    my_classic_theme
  
  return(list(p1, p2, p3)) 
} 

# Create the vpc-plots for model 1
path_to_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_to_vpc_data <- "../../Intermediate/VPC_model1.csv"
plot_list <- create_vpc_plots(path_to_data, path_to_vpc_data)  
p1 <- plot_list[[1]]  
p2 <- plot_list[[2]]  
p3 <- plot_list[[3]]  
ggsave("../../Result/Figures/Model1_nlme/VPC_30_70.pdf", p1, width = 9, height = 6)
ggsave("../../Result/Figures/Model1_nlme/VPC_20_80.pdf", p2, width = 9, height = 6)
ggsave("../../Result/Figures/Model1_nlme/VPC_05_95.pdf", p3, width = 9, height = 6)

# Create the vpc-plots for model 2
path_to_data <- "../../Intermediate/Data_whole_tidy_filt.csv"
path_to_vpc_data <- "../../Intermediate/VPC_model2.csv"
plot_list <- create_vpc_plots(path_to_data, path_to_vpc_data)  
p1 <- plot_list[[1]]  
p2 <- plot_list[[2]]  
p3 <- plot_list[[3]]  
ggsave("../../Result/Figures/Model2_nlme/VPC_30_70_mod2.pdf", p1, width = 9, height = 6)
ggsave("../../Result/Figures/Model2_nlme/VPC_20_80_mod2.pdf", p2, width = 9, height = 6)
ggsave("../../Result/Figures/Model2_nlme/VPC_05_95_mod2.pdf", p3, width = 9, height = 6)