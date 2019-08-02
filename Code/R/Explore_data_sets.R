library(tidyverse)
library(RColorBrewer)
library(latex2exp)

RNGkind("L'Ecuyer-CMRG")

# General plotting parameters 
my_min_theme <- theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13))
my_classic_theme <- theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                            plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=13)) + theme(panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                                                                  colour = "grey"))

# For certain kinds of plots where one parameter increses 
my_palette <- brewer.pal(n = 9, name = "GnBu")[-c(1, 2, 3)]

# colour-blind friendly palette
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The aim with this file is to explore the Suc2 data when the glucose level is decerased from 4 % to 
# 0.5 % at time zero. Overall there are two data set. From an experimental setup they are basically identical, 
# but they were performed at different days. This file will process these data-sets. More specific it will create graphs 
# of the data sets, filter each of them and in the end combine the two data sets.
# Lastly, this script should be run togehter with the run all script in order to have a correct directory structure, else 
# error are likelly to be encountered. 
# It should be noted that the data goes through two rounds of filtering. In the first round the most obvious extreme 
# observatoions are removed. In the second round observations that show very strange behaviour in individual plots 
# are removed. 

# Read the data 
suc2_data1 <- read_csv("../../Data/Suc2_0d1Glc_1.csv") %>% 
  select(-X1) %>% 
  rename("X1" = "X1_1")

suc2_data2 <- read_csv("../../Data/Suc2_0d1Glc_2.csv") %>% 
  select(-X1) %>% 
  rename("X1" = "X1_1")

# -----------------------------------------------------------------------------------------------------------------
# Explore suc2_data1 
# -----------------------------------------------------------------------------------------------------------------
# Filter out the mean columns only and add the timestamps and make data tidy 
data1_mean <- suc2_data1 %>% 
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = dim(suc2_data1)[1])) %>% 
  select(t, everything()) %>% 
  gather(starts_with("Mean"), value = "Mean", key = "Cell") %>% 
  mutate(Cell = as.factor(Cell))

# Calculate the mean and median for each column
mean_median1 <- data1_mean %>% 
  group_by(t) %>% 
  summarise(mean = mean(Mean), 
            median = median(Mean))

# Plot the data 
pdf(file = "../../Result/Figures/Data_set/Suc2_data1_line_plot.pdf", height = 5, width = 9)
cols = c("Mean" = cbPalette[4], "Median" = cbPalette[6])
p1_line <- ggplot(data1_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = mean_median1, aes(t, mean, color = "Mean"), size = 1.5) + 
  geom_line(data = mean_median1, aes(t, median, color = "Median"), size = 1.5) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
  scale_color_manual(name = "Measure", values = cols) + 
  my_min_theme
p1_line
garbage <- dev.off()

# Boxplot yields nice overview 
pdf(file = "../../Result/Figures/Data_set/Suc2_data1_box_plot.pdf", height = 5, width = 9)
p1_box <- ggplot(data1_mean, aes(t, Mean)) + 
  geom_boxplot(aes(group = t)) + 
  labs(x = "Time [min]", y = "Intensity") +
  my_min_theme
p1_box
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Explore suc2_data2 
# -----------------------------------------------------------------------------------------------------------------
# Filter out the mean columns only and add the timestamps and make data tidy 
data2_mean <- suc2_data2 %>% 
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = dim(suc2_data2)[1])) %>% 
  select(t, everything()) %>% 
  gather(starts_with("Mean"), value = "Mean", key = "Cell") %>% 
  mutate(Cell = as.factor(Cell))

# Calculate the mean and median for each column
mean_median2 <- data2_mean %>% 
  group_by(t) %>% 
  summarise(mean = mean(Mean), 
            median = median(Mean))

# Plot the data 
pdf(file = "../../Result/Figures/Data_set/Suc2_data2_line_plot.pdf", height = 5, width = 9)
cols = c("Mean" = cbPalette[4], "Median" = cbPalette[6])
p2_line <- ggplot(data2_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = mean_median2, aes(t, mean, color = "Mean"), size = 1.5) + 
  geom_line(data = mean_median2, aes(t, median, color = "Median"), size = 1.5) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
  scale_color_manual(name = "Measure", values = cols) + 
  my_min_theme
p2_line
garbage <- dev.off()

# Boxplot yields nice overview 
pdf(file = "../../Result/Figures/Data_set/Suc2_data2_box_plot.pdf", height = 5, width = 9)
p2_box <- ggplot(data2_mean, aes(t, Mean)) + 
  geom_boxplot(aes(group = t)) + 
  labs(x = "Time [min]", y = "Intensity") + 
  my_min_theme
p2_box 
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Color most extreme data1 to find outliers 
# -----------------------------------------------------------------------------------------------------------------
# Get the most extreme inital intensities 
data1_initial <- data1_mean %>% filter(t == 0) %>% 
  filter(Mean > 600)
data1_end <- data1_mean %>% filter(t == 480) %>%
  filter(Mean > 1200 | Mean < 380)

# Marking the most extreme cells 
data_extreme1 <- data1_mean %>% 
  filter(Cell == "Mean3" | Cell == "Mean5" | Cell == "Mean29" | Cell == "Mean4" | Cell == "Mean39" | Cell == "Mean26"
         | Cell == "Mean6" | Cell == "Mean9" | Cell == "Mean30" | Cell == "Mean70" | Cell == "Mean31")

# Write result to file 
pdf(file = "../../Result/Figures/Data_set/Suc2_data1_line_plot_ext.pdf", height = 5, width = 9)
p1_ext <- ggplot(data1_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = data_extreme1, aes(t, Mean, color = Cell), size = 1.5) + 
  geom_line(data = mean_median1, aes(t, mean), size = 1.5, linetype = "dashed") + 
  geom_line(data = mean_median1, aes(t, median), size = 1.5) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), title = "Extreme observations data 1") +
  scale_color_brewer(palette = "Paired") + 
  my_min_theme
p1_ext
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Color most extreme data2 to find outliers 
# -----------------------------------------------------------------------------------------------------------------
# Get the most extreme inital intensities 
data2_initial <- data2_mean %>% filter(t == 0) %>% 
  filter(Mean > 800)

# Marking the most extreme cells 
data_extreme2 <- data2_mean %>% 
  filter(Cell == "Mean4" | Cell == "Mean15" | Cell == "Mean16" | Cell == "Mean17" | Cell == "Mean20" | Cell == "Mean25"
         | Cell == "Mean45" | Cell == "Mean46" | Cell == "Mean47" | Cell == "Mean83" | Cell == "Mean88" 
         | Cell == "Mean89") %>%
  mutate(Cell = as.factor(as.character(Cell)))

# Write result to file 
pdf(file = "../../Result/Figures/Data_set/Suc2_data2_line_plot_ext.pdf", height = 5, width = 9)
p2_ext <- ggplot(data2_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = data_extreme2, aes(t, Mean, color = Cell), size = 1.5) + 
  geom_line(data = mean_median2, aes(t, mean), size = 1.5, linetype = "dashed") + 
  geom_line(data = mean_median2, aes(t, median), size = 1.5) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]"), title = "Extreme observations data 2") +
  scale_color_brewer(palette = "Paired") + 
  my_min_theme
p2_ext
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Check sizes of the cells for data 1 to find outliers 
# -----------------------------------------------------------------------------------------------------------------
size_data1 <- suc2_data1 %>% 
  select(starts_with("Area")) %>%
  filter(row_number() == 1) %>%
  gather(starts_with("Area"), value = "Area", key = "Cell") %>%
  mutate(Cell_index = 1:70) %>%
  arrange_at("Area", desc) 

pdf(file = "../../Result/Figures/Data_set/Suc2_data1_size.pdf", height = 5, width = 9)
p1_size <- ggplot(size_data1, aes(Cell_index, Area)) + 
  geom_point() + 
  labs(x = "Cell index", y = "Cell area", title = "Data1 cell sizes") + 
  scale_y_log10() + 
  my_min_theme
p1_size 
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Check sizes of the cells for data 2 to find outliers 
# -----------------------------------------------------------------------------------------------------------------
size_data2 <- suc2_data2 %>% 
  select(starts_with("Area")) %>%
  filter(row_number() == 1) %>%
  gather(starts_with("Area"), value = "Area", key = "Cell") %>%
  mutate(Cell_index = 1:130) %>%
  arrange_at("Area", desc) 

pdf(file = "../../Result/Figures/Data_set/Suc2_data2_size.pdf", height = 5, width = 9)
p2_size <- ggplot(size_data2, aes(Cell_index, Area)) + 
  geom_point() + 
  scale_y_log10() +
  labs(x = "Cell index", y = "Cell area", title = "Data2 cell sizes") + 
  my_min_theme
p2_size
garbage <- dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Plot some of the more extreme cells 
# -----------------------------------------------------------------------------------------------------------------
# Color the more extreme 
data_extreme1 <- data1_mean %>% 
  filter(Cell == "Mean3" | Cell == "Mean5" | Cell == "Mean29" | Cell == "Mean4" | Cell == "Mean39" | Cell == "Mean26"
         | Cell == "Mean6" | Cell == "Mean9" | Cell == "Mean30" | Cell == "Mean70" | Cell == "Mean31")

# Write result to file 
pdf(file = "../../Result/Figures/Data_set/Suc2_data1_line_plot_ext.pdf", height = 5, width = 9)
p1_ext <- ggplot(data1_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = data_extreme1, aes(t, Mean, color = Cell), size = 1.0) + 
  geom_line(data = mean_median1, aes(t, mean), size = 1.1, linetype = "dashed") + 
  geom_line(data = mean_median1, aes(t, median), size = 1.1) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
  scale_color_brewer(palette = "Paired") + 
  my_classic_theme + 
  theme(legend.position = "none")
p1_ext
garbage <- dev.off()

data_extreme2 <- data2_mean %>% 
  filter(Cell == "Mean4" | Cell == "Mean15" | Cell == "Mean16" | Cell == "Mean17" | Cell == "Mean20" | Cell == "Mean25"
         | Cell == "Mean45" | Cell == "Mean46" | Cell == "Mean47" | Cell == "Mean83" | Cell == "Mean88" 
         | Cell == "Mean89") %>%
  mutate(Cell = as.factor(as.character(Cell)))

# Write result to file 
pdf(file = "../../Result/Figures/Data_set/Suc2_data2_line_plot_ext.pdf", height = 5, width = 9)
p2_ext <- ggplot(data2_mean, aes(t, Mean)) + 
  geom_line(aes(group = Cell), size = 0.3, color = cbPalette[1]) + 
  geom_line(data = data_extreme2, aes(t, Mean, color = Cell), size = 1.0) + 
  geom_line(data = mean_median2, aes(t, mean), size = 1.1, linetype = "dashed") + 
  geom_line(data = mean_median2, aes(t, median), size = 1.1) + 
  labs(x = "Time [min]", y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
  scale_color_brewer(palette = "Paired") + 
  my_classic_theme + 
  theme(legend.position = "none")
p2_ext
dev.off()


# -----------------------------------------------------------------------------------------------------------------
# Filter both data sets for the most extreme observations
# -----------------------------------------------------------------------------------------------------------------
# Will use regular expressions to filter away certain observations 
data1_filtered <- suc2_data1 %>%
  select(-matches("^[a-zA-Z]+9$")) %>%
  select(-matches("^[a-zA-Z]+6$")) %>%
  select(-matches("^[a-zA-Z]+30$")) %>%
  select(-matches("^[a-zA-Z]+3$")) %>%
  select(-matches("^[a-zA-Z]+70$")) %>%
  select(-matches("^[a-zA-Z]+31$")) %>%
  select(-matches("^[a-zA-Z]+4$")) %>%
  select(-matches("^[a-zA-Z]+26$")) %>%
  select(-matches("^[a-zA-Z]+29$")) 

data2_filtered <- suc2_data2 %>%
  select(-matches("^[a-zA-Z]+46$")) %>%
  select(-matches("^[a-zA-Z]+45$")) %>%
  select(-matches("^[a-zA-Z]+17$")) %>%
  select(-matches("^[a-zA-Z]+20$")) %>%
  select(-matches("^[a-zA-Z]+88$")) %>%
  select(-matches("^[a-zA-Z]+130$")) %>%
  select(-matches("^[a-zA-Z]+16$")) %>%
  select(-matches("^[a-zA-Z]+15$")) %>%
  select(-matches("^[a-zA-Z]+25$")) %>%
  select(-matches("^[a-zA-Z]+47$")) %>%
  select(-matches("^[a-zA-Z]+83$")) %>%
  select(-matches("^[a-zA-Z]+89$")) 


# -----------------------------------------------------------------------------------------------------------------
# Make the data tidy as this is the format required by monolix 
# -----------------------------------------------------------------------------------------------------------------
## Data 1 
data1_tidy_tmp <- data1_filtered %>%
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = 97)) %>%
  select(t, everything()) %>%
  gather(starts_with("Mean"), value = "Mean", key = "Cell_index") 
data1_tidy_area_tmp <- data1_filtered %>%
  select(starts_with("Area")) %>% 
  gather(starts_with("Area"), value = "Area", key = "Cell_index_area") 
n_cells <- as.integer(dim(data1_tidy_tmp)[1] / 97)
data1_tidy <- data1_tidy_tmp %>%
  bind_cols(data1_tidy_area_tmp) %>%
  select(-"Cell_index_area") %>%
  mutate(Cell = as.factor(rep(1:n_cells, each = 97))) %>%
  mutate(t = t / max(t)) %>%
  mutate(Intensity = Mean / 100) 

## Data 2
data2_tidy_tmp <- data2_filtered %>%
  select(starts_with("Mean")) %>% 
  mutate(t = seq(from = 0, by = 5, length.out = 97)) %>%
  select(t, everything()) %>%
  gather(starts_with("Mean"), value = "Mean", key = "Cell_index") 
data2_tidy_area_tmp <- data2_filtered %>%
  select(starts_with("Area")) %>% 
  gather(starts_with("Area"), value = "Area", key = "Cell_index_area") 
n_cells <- as.integer(dim(data2_tidy_tmp)[1] / 97)
data2_tidy <- data2_tidy_tmp %>%
  bind_cols(data2_tidy_area_tmp) %>%
  select(-"Cell_index_area") %>%
  mutate(Cell = as.factor(rep(1:n_cells, each = 97))) %>%
  mutate(t = t / max(t)) %>%
  mutate(Intensity = Mean / 100) 

# Second round of filtering after looking at each individual cell. See report for more details about filtered cells 
# Filter data1 and remove strange cells
data1_tidy_filtered <- data1_tidy %>%
  filter(Cell != 7) %>%
  filter(Cell != 16) %>%
  filter(Cell != 21) %>%
  filter(Cell != 38) %>%
  filter(Cell != 44) %>%
  filter(Cell != 49) %>%
  filter(Cell != 54) %>%
  filter(Cell != 59) 

  # Filter the second data set 
data2_tidy_filtered <- data2_tidy %>%
  filter(Cell != 3) %>%
  filter(Cell != 4)%>%
  filter(Cell != 5) %>%
  filter(Cell != 7) %>%
  filter(Cell != 11) %>%
  filter(Cell != 12) %>%
  filter(Cell != 15) %>%
  filter(Cell != 16) %>%
  filter(Cell != 17) %>%
  filter(Cell != 18) %>%
  filter(Cell != 24) %>%
  filter(Cell != 28) %>%
  filter(Cell != 29) %>%
  filter(Cell != 30) %>%
  filter(Cell != 31) %>%
  filter(Cell != 32) %>%
  filter(Cell != 33) %>%
  filter(Cell != 39) %>%
  filter(Cell != 40) %>%
  filter(Cell != 41) %>%
  filter(Cell != 48) %>%
  filter(Cell != 53) %>%
  filter(Cell != 56) %>%
  filter(Cell != 60) %>%
  filter(Cell != 65) %>%
  filter(Cell != 73) %>%
  filter(Cell != 74) %>%
  filter(Cell != 75) %>%
  filter(Cell != 76) %>%
  filter(Cell != 77) %>%
  filter(Cell != 78) %>%
  filter(Cell != 79) %>%
  filter(Cell != 81) %>%
  filter(Cell != 82) %>%
  filter(Cell != 84) %>%
  filter(Cell != 86) %>%
  filter(Cell != 90) %>%
  filter(Cell != 95) %>%
  filter(Cell != 98) %>%
  filter(Cell != 101) %>%
  filter(Cell != 108) %>%
  filter(Cell != 110) %>%
  filter(Cell != 113) %>%
  filter(Cell != 114) %>%
  filter(Cell != 116) %>%
  filter(Cell != 117) 

# -----------------------------------------------------------------------------------------------------------------
# Create the whole data set in tidy format 
# -----------------------------------------------------------------------------------------------------------------
n_cells1 <- dim(data1_tidy_filtered)[1] / 97
n_cells2 <- dim(data2_tidy_filtered)[1] / 97
n_obs1 <- dim(data1_tidy_filtered)[1] 
n_obs2 <- dim(data2_tidy_filtered)[1]
whole_data_filt_tidy <- data1_tidy_filtered %>%
  bind_rows(data2_tidy_filtered) %>%
  mutate(Experiment = as.factor(c(rep(1, n_obs1), rep(2, n_obs2)))) %>%
  mutate(ID = as.factor(rep(1:(n_cells1 + n_cells2), each = 97))) %>%
  rename("Suc2" = "Intensity")

# Plot the whole data set in a line plot 
mean_median <- whole_data_filt_tidy %>%
  group_by(t) %>%
  summarise(mean = mean(Suc2), 
            median = median(Suc2))
cols = c("Mean" = cbPalette[4], "Median" = cbPalette[6])
p_line_whole <- ggplot(whole_data_filt_tidy, aes(t, Suc2)) + 
  geom_line(aes(group = ID), color = cbPalette[1], size = 0.3) + 
  geom_line(data = mean_median, aes(t, mean, color = "Mean"), size = 1.5) + 
  geom_line(data = mean_median, aes(t, median, color = "Median"), size = 1.5) + 
  labs(x = TeX("Scaled time"), y = TeX("Suc2 intensity \\[A.U.$\\times 10^{-2}$\\]")) +
  scale_color_manual(name = "Measure", values = cols) + 
  my_classic_theme

p_line_whole
ggsave("../../Result/Figures/Data_set/Line_whole_data_set_line.pdf", width = 9, height = 5)

# Plot the different data sets 
p_line_exp <- ggplot(whole_data_filt_tidy, aes(t, Suc2, color = Experiment)) + 
  geom_line(aes(group = ID), size = 0.3) + 
  scale_color_manual(values = cbPalette[-1]) + 
  my_classic_theme
p_line_exp
ggsave("../../Result/Figures/Data_set/Line_whole_data_set_exp_col.pdf", width = 9, height = 5)

# Box plot of both data sets
p_both_box <- ggplot(whole_data_filt_tidy, aes(as.factor(t), Suc2, fill = Experiment)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette[-1]) +
  my_classic_theme +
  theme(axis.text.x = element_blank())
p_both_box
ggsave("../../Result/Figures/Data_set/Box_whole_data_set_exp_col.pdf", width = 9, height = 5)

# Write the whole tidy data set to file 
if(!file.exists("../../Intermediate/Data_whole_tidy_filt.csv")){
  write_csv(whole_data_filt_tidy, "../../Intermediate/Data_whole_tidy_filt.csv")
  Sys.chmod("../../Intermediate/Data_whole_tidy_filt.csv", mode = "444")
}

