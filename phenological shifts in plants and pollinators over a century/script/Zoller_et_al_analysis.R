#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
# This script contains code for reproducing the analyses presented in the Manuscript:

# Zoller L., Vázquez D.P. & Resasco J. Phenological shifts in plants and pollinators over a 
# century disrupt interaction persistence
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________


## Set up ----

# add packages
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggridges)
library(ggplot2)
library(reshape2)
library(bipartite)
library(grid)
library(gridGraphics)
library(cowplot)
library(pracma)
library(performance)

# load data
# Pollinators
Poll_data <- read.csv("data/Data_eleSubset_PollSubset_R.csv")
str(Poll_data)

# Repeat rows based on the 'Freq' variable
Poll_data <- Poll_data[rep(row.names(Poll_data), times = Poll_data$Freq), ]
Poll_data$Freq <- 1

# Reorder Species by mean Julian_day
Poll_data$Poll_spec <- reorder(Poll_data$Poll_spec,
                               ave(Poll_data$Julian_day,
                                   Poll_data$Poll_spec,
                                   FUN = mean))
head(Poll_data)


# Plants
Plant_data <- read.csv("data/Data_eleSubset_PlantSubset_R.csv")
str(Plant_data)

# Repeat rows based on the 'Freq' variable
Plant_data <- Plant_data[rep(row.names(Plant_data), times = Plant_data$Freq), ]
Plant_data$Freq <- 1

# Reorder Species by mean Julian_day
Plant_data$Plant_spec <- reorder(Plant_data$Plant_spec,
                                 ave(Plant_data$Julian_day,
                                     Plant_data$Plant_spec,
                                     FUN = mean))
head(Plant_data)

#_______________________________________________________________________

# Part 1: Shifts in phenology - plant and pollinator community ----
#_______________________________________________________________________

### Calculate flowering and flight metrics ----

# define function to calculate onset, peak, end and duration of flight/flowering
calculate_metrics <- function(data) {
  metrics <- data %>%
    summarize(
      Peak_julian_day = Julian_day[which.max(Freq)],
      Onset_julian_day = min(Julian_day),
      End_julian_day = max(Julian_day),
      Duration = max(Julian_day) - min(Julian_day),
      .groups = 'drop'
    )
  return(metrics)
}

Plant_metrics <- Plant_data %>%
  group_by(Plant_spec, Period) %>%
  do(calculate_metrics(.))

Poll_metrics <- Poll_data %>%
  group_by(Poll_spec, Period) %>%
  do(calculate_metrics(.))


# Add a Type column to each dataframe
Plant_metrics <- Plant_metrics %>% 
  mutate(Type = "Plants")

Poll_metrics <- Poll_metrics %>% 
  mutate(Type = "Pollinators")

# Combine the two dataframes
combined_metrics <- bind_rows(Plant_metrics, Poll_metrics)

# Calculate mean and standard error separately for Past and Present periods
mean_se_period <- combined_metrics %>%
  group_by(Type, Period) %>%
  summarize(
    mean_onset = mean(Onset_julian_day),
    mean_peak = mean(Peak_julian_day),
    mean_end = mean(End_julian_day),
    mean_duration = mean(Duration),
    se_onset = sd(Onset_julian_day) / sqrt(n()),
    se_peak = sd(Peak_julian_day) / sqrt(n()),
    se_end = sd(End_julian_day) / sqrt(n()),
    se_duration = sd(Duration) / sqrt(n()),
    .groups = 'drop'  # to ungroup after summarizing
  )

mean_se_period


### Perform t-tests and visualize ----

#### Onset ----

A <- ggplot(mean_se_period, aes(x = Type, y = mean_onset, color = Period, group = Period)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_onset - se_onset, ymax = mean_onset + se_onset), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Past" = "#66787A", "Present" = "skyblue")) +
  labs(x = "",
       y = "Mean onset (day of year)") +
  theme_minimal() +
  ylim(140, 250) +
  theme(legend.position = "none") 
A

# T-tests
past_data <- combined_metrics %>% filter(Period == "Past")
present_data <- combined_metrics %>% filter(Period == "Present")
poll_data <- combined_metrics %>% filter(Type == "Pollinators")
plant_data <- combined_metrics %>% filter(Type == "Plants")

(t_test_past <- t.test(Onset_julian_day ~ Type, data = past_data))
(t_test_poll <- t.test(Onset_julian_day ~ Period, data = poll_data))
(t_test_present <- t.test(Onset_julian_day ~ Type, data = present_data))
(t_test_plant <- t.test(Onset_julian_day ~ Period, data = plant_data))

#### Peak ----

B <- ggplot(mean_se_period, aes(x = Type, y = mean_peak, color = Period, group = Period)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_peak - se_peak, ymax = mean_peak + se_peak), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Past" = "#66787A", "Present" = "skyblue")) +
  labs(x = "",
       y = "Mean peak (day of year)") +
  theme_minimal() +
  ylim(140, 250) +
  theme(legend.position = "none") 
B

# T-tests
(t_test_past <- t.test(Peak_julian_day ~ Type, data = past_data))
(t_test_poll <- t.test(Peak_julian_day ~ Period, data = poll_data))
(t_test_present <- t.test(Peak_julian_day ~ Type, data = present_data))
(t_test_plant <- t.test(Peak_julian_day ~ Period, data = plant_data))

#### End ----

C <- ggplot(mean_se_period, aes(x = Type, y = mean_end, color = Period, group = Period)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_end - se_end, ymax = mean_end + se_end), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Past" = "#66787A", "Present" = "skyblue")) +
  labs(x = "",
       y = "Mean end (day of year)") +
  theme_minimal()+
  ylim(140, 250) +
  theme(legend.position = "none") 
C

# T-tests
(t_test_past <- t.test(End_julian_day ~ Type, data = past_data))
(t_test_poll <- t.test(End_julian_day ~ Period, data = poll_data))
(t_test_present <- t.test(End_julian_day ~ Type, data = present_data))
(t_test_plant <- t.test(End_julian_day ~ Period, data = plant_data))


#### Duration ----

D <- ggplot(mean_se_period, aes(x = Type, y = mean_duration, color = Period, group = Period)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_duration - se_duration, ymax = mean_duration + se_duration), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Past" = "#66787A", "Present" = "skyblue")) +
  labs(x = "",
       y = "Mean duration (days)") +
  theme_minimal() +
  ylim(20, 65)
D

# T-tests
(t_test_past <- t.test(Duration ~ Type, data = past_data))
(t_test_poll <- t.test(Duration ~ Period, data = poll_data))
(t_test_present <- t.test(Duration ~ Type, data = present_data))
(t_test_plant <- t.test(Duration ~ Period, data = plant_data))

# arrange plots
cowplot::plot_grid(A, B, C, D, 
                   labels = "AUTO",
                   ncol = 2, nrow = 2,
                   label_size = 12)





#_______________________________________________________________________

# Part 2: Shifts in pollinator phenology - individual species ----
#_______________________________________________________________________

### Visualize ridge plot  of pollinator flight activity ----

# Reorder Species by median Julian_day
Poll_data$Poll_spec <- reorder(Poll_data$Poll_spec,
                               ave(Poll_data$Julian_day,
                                   Poll_data$Poll_spec,
                                   FUN = mean))

# plot
Ridge_plot_poll <- ggplot(Poll_data, aes(x = Julian_day, y = Poll_spec, fill = Period)) + 
                   geom_density_ridges(scale = 1, 
                                       jittered_points = FALSE, 
                                       point_fill = "grey80", 
                                       alpha = .6,
                                       quantile_lines = FALSE, quantiles = 2) +
                   scale_fill_manual(values = c("#66787A", "#D9af27")) +  
                   labs(x = "Day of year", 
                        y = "Pollinator species") +
                   theme_minimal() +
                   theme(panel.grid.major.x = element_blank(),
                         panel.grid.minor.x = element_blank(),
                         legend.position = "none") 

Ridge_plot_poll


### Visualize point range plot ----

# Calculate mean, standard error, and 95% CI of Julian_day for each species and period
summary_data <- Poll_data %>%
  group_by(Poll_spec, Period) %>%
  summarize(
    mean_julian_day = mean(Julian_day),
    peak_julian_day = Julian_day[which.max(Freq)],
    sd_julian_day = sd(Julian_day),
    se_julian_day = sd_julian_day / sqrt(n()),
    ci_lower_mean = mean_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_mean = mean_julian_day + qt(0.975, df = n() - 1) * se_julian_day,
    ci_lower_peak = peak_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_peak = peak_julian_day + qt(0.975, df = n() - 1) * se_julian_day,
    first_observation = min(Julian_day),
    last_observation = max(Julian_day)
  )
head(summary_data)

Point_plot_poll <- ggplot(summary_data, aes(x = mean_julian_day, 
                                            y = Poll_spec, 
                                            color = Period, 
                                            group = Period)) +
  geom_pointrange(aes(xmin = ci_lower_mean, xmax = ci_upper_mean),
                  position = position_dodge(0.5)) +
  scale_color_manual(values = c("#66787A", "#D9af27")) + 
  theme_minimal() +
  coord_cartesian(xlim = c(120, 260)) +
  labs(x = "Day of year") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none", 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank()) 

Point_plot_poll


### Perform T-tests for shift in mean day of year ----

species_list <- unique(Poll_data$Poll_spec)

# Create an empty data frame to store results
poll_results <- data.frame(Species = character(),
                      poll_MFD_past = numeric(),
                      poll_MFD_pres = numeric(),
                      poll_MFD_diff = numeric(),
                      poll_MinFD_past = numeric(),
                      poll_MinFD_pres = numeric(),
                      poll_MaxFD_past = numeric(),
                      poll_MaxFD_pres = numeric(),
                      poll_Test_Statistic = numeric(),
                      poll_DegFree = numeric(),
                      poll_p_value = numeric(),
                      poll_Stderr = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each species
for (species in species_list) {
  # Subset data for the current species
  species_data_past <- subset(Poll_data, Period == "Past" & Poll_spec == species)$Julian_day
  species_data_present <- subset(Poll_data, Period == "Present" & Poll_spec == species)$Julian_day
  
  # Calculate mean Julian day for each period
  mean_julian_day_past <- mean(species_data_past)
  mean_julian_day_present <- mean(species_data_present)
  mean_julian_day_difference <- mean_julian_day_present - mean_julian_day_past
  
  # Calculate earliest Julian day for each period
  min_julian_day_past <- min(species_data_past)
  min_julian_day_present <- min(species_data_present)
  
  # Calculate latest Julian day for each period
  max_julian_day_past <- max(species_data_past)
  max_julian_day_present <- max(species_data_present)
  
  # Perform t-test
  t_test_species <- t.test(species_data_past, species_data_present)
  
  # Store test results
  poll_results <- rbind(poll_results, data.frame(Species = species,
                                       poll_MFD_past = mean_julian_day_past,
                                       poll_MFD_pres = mean_julian_day_present,
                                       poll_MinFD_past = min_julian_day_past,
                                       poll_MinFD_pres = min_julian_day_present,
                                       poll_MFD_diff = mean_julian_day_difference,
                                       poll_MaxFD_past = max_julian_day_past,
                                       poll_MaxFD_pres = max_julian_day_present,
                                       poll_Test_Statistic = t_test_species$statistic,
                                       poll_DegFree = t_test_species$parameter,
                                       poll_p_value = t_test_species$p.value,
                                       poll_Stderr = t_test_species$stderr,
                                       stringsAsFactors = FALSE))
}

# write.csv(poll_results, file = "data/poll_results.csv", row.names = FALSE)

# define custom plotting order
custom_order <- c("Colletes paniscus", "Andrena prunorum", 
                  "Osmia pentstemonis", "Pseudomasaris vespoides",
                  "Andrena crataegi", "Osmia densa",
                  "Bombus huntii","Hylaeus basalis",
                  "Anthophora ursina", "Hoplitis albifrons",
                  "Apis mellifera", "Megachile melanophaea", 
                  "Osmia proxima", "Osmia bruneri",
                  "Bombus occidentalis", "Eristalis stipator",
                  "Anthophora terminalis", "Hylaeus annulatus",
                  "Megachile relativa", "Bombus centralis",
                  "Hylaeus wootoni", "Bombus bifarius",
                  "Bombus appositus", "Poanes taxiles", 
                  "Bombus melanopygus")

# Convert Species column to factor with custom order
poll_results$Species <- factor(poll_results$Species, levels = custom_order)

# plot
T_test_poll <- ggplot(poll_results, aes(x = poll_MFD_pres - poll_MFD_past, y = Species)) +
  geom_point(aes(color = ifelse(poll_p_value        < 0.05, "Significant", "Non-significant")), size = 2) +  
  geom_errorbar(aes(xmin = poll_MFD_pres - poll_MFD_past - 1.96 * poll_Stderr,
                    xmax = poll_MFD_pres - poll_MFD_past + 1.96 * poll_Stderr),
                width = 0.2, color = "black", alpha = 0.5) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Δ mean Julian day", y = "") +  
  scale_color_manual(name = "Significance",
                     values = c("Significant" = "red", "Non-significant" = "black"),
                     labels = c("Non-significant", "Significant")) + 
  guides(color = guide_legend(title = "Significance")) + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.position = "none", 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank()) #+

T_test_poll


# Arrange plots 
cowplot::plot_grid(
  plot_grid(Ridge_plot_poll, ncol = 1, rel_widths = c(5, 5)),  
  plot_grid(Point_plot_poll, T_test_poll, ncol = 2, rel_widths = c(0.1, 0.1)), 
  align = "hv",
  ncol = 2  ,
  label_size = 12,
  labels = "AUTO"
)





#_______________________________________________________________________

# Part 3: Shifts in plant phenology - individual species ----
#_______________________________________________________________________

### Visualize as ridge plot ----

# Reorder Species by median Julian_day
Plant_data$plant_spec <- reorder(Plant_data$Plant_spec,
                                 ave(Plant_data$Julian_day,
                                     Plant_data$Plant_spec,
                                     FUN = mean))


# plot
Ridge_plot_plant <- ggplot(Plant_data, aes(x = Julian_day, y = plant_spec, fill = Period)) + 
  geom_density_ridges(scale = 1, 
                      jittered_points = FALSE, 
                      point_fill = "grey80", 
                      alpha = .6,
                      quantile_lines = FALSE, quantiles = 2) +
  scale_fill_manual(values = c("#66787A", "#00CC00")) +  
  labs(x = "Day of year", 
       y = "Pollinator species") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none") 

Ridge_plot_plant


### Visualize as point range plot ----

# Calculate mean, standard error, and 95% CI of Julian_day for each Species and Period
summary_data <- Plant_data %>%
  group_by(Plant_spec, Period) %>%
  summarize(
    mean_julian_day = mean(Julian_day),
    peak_julian_day = Julian_day[which.max(Freq)],
    sd_julian_day = sd(Julian_day),
    se_julian_day = sd_julian_day / sqrt(n()),
    ci_lower_mean = mean_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_mean = mean_julian_day + qt(0.975, df = n() - 1) * se_julian_day,
    ci_lower_peak = peak_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_peak = peak_julian_day + qt(0.975, df = n() - 1) * se_julian_day,
    first_observation = min(Julian_day),
    last_observation = max(Julian_day)
  )

head(summary_data)

custom_order <- c("Rubus deliciosus",
                  "Penstemon virens", "Penstemon secundiflorus",
                  "Rubus idaeus var. strigosus", "Penstemon glaber",
                  "Geranium richardsonii", "Aquilegia coerulea",
                  "Geranium caespitosum", "Aconitum columbianum", 
                  "Chamaenerion angustifolium")


# Convert Species column to factor with custom order
summary_data$Plant_spec <- factor(summary_data$Plant_spec, levels = custom_order)

Point_plot_plant <- ggplot(summary_data, aes(x = mean_julian_day, 
                                             y = Plant_spec, 
                                             color = Period, 
                                             group = Period)) +
  geom_pointrange(aes(xmin = ci_lower_mean, xmax = ci_upper_mean),
                  position = position_dodge(0.5)) +
  scale_color_manual(values = c("#66787A", "#00CC00")) + 
  theme_minimal() +
  coord_cartesian(xlim = c(120, 260)) +
  labs(x = "Day of year") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none", 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank()
  ) 

Point_plot_plant


### Perform T-tests for shift in mean day of year ----

species_list <- unique(Plant_data$Plant_spec)

# Create an empty data frame to store results
plant_results <- data.frame(Species = character(),
                            plant_MFD_past = numeric(),
                            plant_MFD_pres = numeric(),
                            plant_MFD_diff = numeric(),
                            plant_MinFD_past = numeric(),
                            plant_MinFD_pres = numeric(),
                            plant_MaxFD_past = numeric(),
                            plant_MaxFD_pres = numeric(),
                            plant_Test_Statistic = numeric(),
                            plant_DegFree = numeric(),
                            plant_p_value = numeric(),
                            plant_Stderr = numeric(),
                            stringsAsFactors = FALSE)

# Loop through each species
for (species in species_list) {
  # Subset data for the current species
  species_data_past <- subset(Plant_data, Period == "Past" & Plant_spec == species)$Julian_day
  species_data_present <- subset(Plant_data, Period == "Present" & Plant_spec == species)$Julian_day
  
  # Calculate mean Julian day for each period
  mean_julian_day_past <- mean(species_data_past)
  mean_julian_day_present <- mean(species_data_present)
  mean_julian_day_difference <- mean_julian_day_present - mean_julian_day_past
  
  # Calculate earliest Julian day for each period
  min_julian_day_past <- min(species_data_past)
  min_julian_day_present <- min(species_data_present)
  
  # Calculate latest Julian day for each period
  max_julian_day_past <- max(species_data_past)
  max_julian_day_present <- max(species_data_present)
  
  # Perform t-test
  t_test_species <- t.test(species_data_past, species_data_present)
  
  # Store test results
  plant_results <- rbind(plant_results, data.frame(Species = species,
                                                   plant_MFD_past = mean_julian_day_past,
                                                   plant_MFD_pres = mean_julian_day_present,
                                                   plant_MinFD_past = min_julian_day_past,
                                                   plant_MinFD_pres = min_julian_day_present,
                                                   plant_MFD_diff = mean_julian_day_difference,
                                                   plant_MaxFD_past = max_julian_day_past,
                                                   plant_MaxFD_pres = max_julian_day_present,
                                                   plant_Test_Statistic = t_test_species$statistic,
                                                   plant_DegFree = t_test_species$parameter,
                                                   plant_p_value = t_test_species$p.value,
                                                   plant_Stderr = t_test_species$stderr,
                                                   stringsAsFactors = FALSE))
}

#write.csv(plant_results, file = "data/plant_results.csv", row.names = FALSE)


custom_order <- c("Rubus deliciosus",
                  "Penstemon virens", "Penstemon secundiflorus",
                  "Rubus idaeus var. strigosus", "Penstemon glaber",
                  "Geranium richardsonii", "Aquilegia coerulea",
                  "Geranium caespitosum", "Aconitum columbianum", 
                  "Chamaenerion angustifolium")


# Convert Species column to factor with custom order
plant_results$Species <- factor(plant_results$Species, levels = custom_order)


# plot
T_test_plant <- ggplot(plant_results, aes(x = plant_MFD_pres - plant_MFD_past, y = Species)) +
  geom_point(aes(color = ifelse(plant_p_value        < 0.05, "Significant", "Non-significant")), size = 2) +  # All points
  geom_errorbar(aes(xmin = plant_MFD_pres - plant_MFD_past - 1.96 * plant_Stderr,
                    xmax = plant_MFD_pres - plant_MFD_past + 1.96 * plant_Stderr),
                width = 0.2, color = "black", alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Δ mean Julian day", y = "") +  
  scale_color_manual(name = "Significance",
                     values = c("Significant" = "red", "Non-significant" = "black"),
                     labels = c("Non-significant", "Significant")) +  
  guides(color = guide_legend(title = "Significance")) +  
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank()) 


T_test_plant

# Arrange plots 
cowplot::plot_grid(
  plot_grid(Ridge_plot_plant, ncol = 1, rel_widths = c(5, 5)),  
  plot_grid(Point_plot_plant, T_test_plant, ncol = 2, rel_widths = c(0.1, 0.1)), 
  align = "hv",
  ncol = 2  ,
  label_size = 12,
  labels = "AUTO"
)





#_______________________________________________________________________

# Part 4: Persistence of interactions ----
#_______________________________________________________________________

# make new variable "Interaction", concatenating Plant and Pollinator
Poll_data <- Poll_data %>% 
  unite("Interaction", Poll_spec, Plant_spec, sep = "-", remove = FALSE)

# calculate interaction persistence
int_persistence <- Poll_data %>%
  group_by(Interaction, Period) %>%
  summarise(Total_Freq = sum(Freq)) %>%
  ungroup() %>%
  complete(Interaction, Period, fill = list(Total_Freq = 0)) %>%
  pivot_wider(names_from = Period, values_from = Total_Freq, names_prefix = "freq_") %>%
  mutate(
    interaction_persistence = ifelse(freq_Past > 0 & freq_Present > 0, 1, 0),
    interaction_gain_loss = freq_Present - freq_Past,
    Poll_species = sub("-.*", "", Interaction),
    Plant_species = sub(".*-", "", Interaction)
  )

head(int_persistence)
str(int_persistence)

# Add information on phenology of plants and pollinators to the data
#poll_results <- read.csv("data/poll_results.csv")
poll_results <- rename(poll_results, Poll_species = Species)
poll_results <- poll_results %>% select(- poll_Test_Statistic, 
                                        - poll_DegFree, 
                                        - poll_p_value, 
                                        - poll_Stderr)
head(poll_results)


#plant_results <- read.csv("data/plant_results.csv")
plant_results <- rename(plant_results, Plant_species = Species)
plant_results <- plant_results %>% select(- plant_Test_Statistic, 
                                          - plant_DegFree, 
                                          - plant_p_value, 
                                          - plant_Stderr)
head(plant_results)

# Combine data frames
data <- merge(int_persistence, plant_results, by = "Plant_species")
data <- merge(data, poll_results, by = "Poll_species")
str(data)

# calculate duration of activity period of plants and pollinators
data$poll_activity_period_past <- data$poll_MaxFD_past - data$poll_MinFD_past
data$poll_activity_period_pres <- data$poll_MaxFD_pres - data$poll_MinFD_pres
data$plant_activity_period_past <- data$plant_MaxFD_past - data$plant_MinFD_past
data$plant_activity_period_pres <- data$plant_MaxFD_pres - data$plant_MinFD_pres

# calculate days of overlap of plants and polliantors:
data$plant_poll_overlap_days_past <- pmax(0, pmin(data$plant_MaxFD_past, data$poll_MaxFD_past) - pmax(data$plant_MinFD_past, data$poll_MinFD_past))
data$plant_poll_overlap_days_present <- pmax(0, pmin(data$plant_MaxFD_pres, data$poll_MaxFD_pres) - pmax(data$plant_MinFD_pres, data$poll_MinFD_pres))
data$overlap_change <- data$plant_poll_overlap_days_present - data$plant_poll_overlap_days_past

# remove rows where interaction were not observed in the past
filtered_data <- data %>% filter(freq_Past > 0)

# Predict probabilities for the logistic regression curve
model <- glm(interaction_persistence ~ overlap_change, data = filtered_data, family = binomial)
summary(model)

pred_data <- data.frame(overlap_change = seq(min(filtered_data$overlap_change), max(filtered_data$overlap_change), length.out = 100))
pred_data$predicted_persistence <- predict(model, newdata = pred_data, type = "response")

# Create the scatterplot with ggplot
ggplot(filtered_data, aes(x = overlap_change, y = interaction_persistence)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, color = "red") +
  labs(x = "Change in days of overlap", y = "interaction persistence") +
  xlim(-40, 70) +
  theme_minimal()

