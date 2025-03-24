#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
# This script contains code for reproducing the analyses and figures presented in the Manuscript:

# Zoller L., Vázquez D.P. & Resasco J. Phenological shifts in plants and pollinators over a 
# century disrupt interaction persistence

# The manuscript is currently under consideration for publication in the American Naturalist
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________


# Set up ----

# add packages
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggridges)
library(ggplot2)
library(bipartite)
library(grid)
library(gridGraphics)
library(cowplot)
library(pracma)
library(performance)
library(lme4)
library(fitdistrplus)
library(boot)

# Load data
# Pollinators
Poll_data <- read.csv("data/Data_eleSubset_PollSubset_R.csv")
str(Poll_data)

# Repeat rows based on the 'Freq' variable
Poll_data <- Poll_data %>%
  slice(rep(row_number(), Freq)) %>%
  mutate(Freq = 1)

head(Poll_data)


# Plants
Plant_data <- read.csv("data/Data_eleSubset_PlantSubset_R.csv")
str(Plant_data)

# Repeat rows based on the 'Freq' variable
Plant_data <- Plant_data %>%
  slice(rep(row_number(), Freq)) %>%
  mutate(Freq = 1)

head(Plant_data)

#_______________________________________________________________________
#_______________________________________________________________________
# Results and Figures presented in main manuscript ----
#_______________________________________________________________________
#_______________________________________________________________________
#_______________________________________________________________________
# 1: Shifts in phenology - plant and pollinator community ----
#_______________________________________________________________________

# Create function to calculate phenological metrics using MinMax and percentile approaches

calculate_metrics <- function(data) {
  
  # MinMax approach
  MinMax_metrics <- data %>%
    summarize(
      Onset = min(Julian_day),
      Mean = mean(Julian_day),
      End = max(Julian_day),
      Duration = max(Julian_day) - min(Julian_day),
      n = n(), 
      .groups = 'drop'
    )
  
  # Percentile-based approach
  dens <- density(data$Julian_day, adjust = 1)
  cdf <- cumsum(dens$y) / sum(dens$y) 
  
  density_quantiles <- tibble(
    Onset_percentile = approx(cdf, dens$x, xout = 0.05)$y,
    Median_percentile = approx(cdf, dens$x, xout = 0.5)$y,
    End_percentile = approx(cdf, dens$x, xout = 0.95)$y,
    Duration_percentile = approx(cdf, dens$x, xout = 0.95)$y - approx(cdf, dens$x, xout = 0.05)$y,
    
  )
  
  # Combine empirical and quantile-based metrics
  bind_cols(MinMax_metrics, density_quantiles)
}

# Apply function to plant and pollinator data
Plant_metrics <- Plant_data %>%
  group_by(Species = Plant_spec, Period) %>%
  group_modify(~ calculate_metrics(.x)) %>%
  ungroup() %>%
  mutate(Type = "Plants")

Poll_metrics <- Poll_data %>%
  group_by(Species = Poll_spec, Period) %>%
  group_modify(~ calculate_metrics(.x)) %>%
  ungroup() %>%
  mutate(Type = "Pollinators")

# Combine plant and pollinator datasets
Metrics <- bind_rows(Plant_metrics, Poll_metrics)

# Compute mean and standard error by type and period
metrics_mean_se <- Metrics %>%
  group_by(Type, Period) %>%
  summarize(
    
    # empirical metrics
    mean_onset = mean(Onset),
    se_onset = sd(Onset) / sqrt(n()),
    mean_mean = mean(Mean),
    se_mean = sd(Mean) / sqrt(n()),
    mean_end_day = mean(End),
    se_end_day = sd(End) / sqrt(n()),
    mean_duration = mean(Duration),
    se_duration = sd(Duration) / sqrt(n()),
    
    # Percentile-based metrics
    mean_onset_perc = mean(Onset_percentile),
    se_onset_perc = sd(Onset_percentile) / sqrt(n()),
    mean_median_perc = mean(Median_percentile),
    se_median_perc = sd(Median_percentile) / sqrt(n()),
    mean_end_perc = mean(End_percentile),
    se_end_perc = sd(End_percentile) / sqrt(n()),
    mean_Duration_percentile = mean(Duration_percentile),           
    se_Duration_percentile = sd(Duration_percentile) / sqrt(n()),
    
    .groups = 'drop'
  )

metrics_mean_se$n <- c(10, 10, 25, 25) # sample size of plants (n=10) and pollinators (n=25)

## Fig. 1 ----
# Point plots showing the means and standard deviations of onset, mean, end and 
# duration of pollinator flight and plant flowering

### (A) Onset ----
Fig1A <- ggplot(metrics_mean_se, aes(x = Period, y = mean_onset, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_onset - se_onset, ymax = mean_onset + se_onset), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "Onset") +
  theme_minimal() +
  ylim(140, 250) 

### (B) Mean ----
Fig1B <- ggplot(metrics_mean_se, aes(x = Period, y = mean_mean, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_mean - se_mean, ymax = mean_mean + se_mean), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "Mean") +
  theme_minimal() +
  ylim(140, 250) 

### (C) End ----
Fig1C <- ggplot(metrics_mean_se, aes(x = Period, y = mean_end_day, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_end_day - se_end_day, ymax = mean_end_day + se_end_day), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "End") +
  theme_minimal()+
  ylim(140, 250) 

### (D) Duration ----
Fig1D <- ggplot(metrics_mean_se, aes(x = Period, y = mean_duration, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_duration - se_duration, ymax = mean_duration + se_duration), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "Duration") +
  theme_minimal() +
  theme(legend.position = "right") +
  ylim(19, 65)

# Arrange plots
Fig1A <- Fig1A + theme(legend.position = "none")
Fig1B <- Fig1B + theme(legend.position = "none")
Fig1C <- Fig1C + theme(legend.position = "none")
Fig1D <- Fig1D + theme(legend.position = "none")

Fig1 <- cowplot::plot_grid(Fig1A, Fig1B, Fig1C, Fig1D, 
                           labels = "AUTO",
                           ncol = 2, nrow = 2,
                           label_size = 12)
#ggsave("Fig_1.pdf", width = 5, height = 6, dpi = 300)


## Statistical tests ----
past_data <- Metrics %>% filter(Period == "Past")
present_data <- Metrics %>% filter(Period == "Present")
poll_data <- Metrics %>% filter(Type == "Pollinators")
plant_data <- Metrics %>% filter(Type == "Plants")

{
  MinMax_test_results <- list()
  # onset
  (MinMax_test_results$t_test_onset_past <- t.test(Onset ~ Type, data = past_data))
  (MinMax_test_results$t_test_onset_poll <- t.test(Onset ~ Period, data = poll_data))
  (MinMax_test_results$t_test_onset_present <- t.test(Onset ~ Type, data = present_data))
  (MinMax_test_results$t_test_onset_plant <- t.test(Onset ~ Period, data = plant_data))
  # mean
  (MinMax_test_results$t_test_mean_past <- t.test(Mean ~ Type, data = past_data))
  (MinMax_test_results$t_tes_mean_poll <- t.test(Mean ~ Period, data = poll_data))
  (MinMax_test_results$t_test_mean_present <- t.test(Mean ~ Type, data = present_data))
  (MinMax_test_results$t_test_mean_plant <- t.test(Mean ~ Period, data = plant_data))
  # end
  (MinMax_test_results$t_test_end_past <- t.test(End ~ Type, data = past_data))
  (MinMax_test_results$t_test_end_poll <- t.test(End ~ Period, data = poll_data))
  (MinMax_test_results$t_test_end_present <- t.test(End ~ Type, data = present_data))
  (MinMax_test_results$t_test_end_plant <- t.test(End ~ Period, data = plant_data))
  # duration
  (MinMax_test_results$t_test_duration_past <- t.test(Duration ~ Type, data = past_data))
  (MinMax_test_results$t_test_duration_poll <- t.test(Duration ~ Period, data = poll_data))
  (MinMax_test_results$t_test_duration_present <- t.test(Duration ~ Type, data = present_data))
  (MinMax_test_results$t_test_duration_plant <- t.test(Duration ~ Period, data = plant_data))
  
  MinMax_test_results
}


#_______________________________________________________________________

# 2: Shifts in pollinator phenology - individual species ----
#_______________________________________________________________________

## Fig. 2 ----
### (A) Ridge plot ----

# Reorder Species by mean Julian_day
Poll_data$Poll_spec <- reorder(Poll_data$Poll_spec,
                               ave(Poll_data$Julian_day,
                                   Poll_data$Poll_spec,
                                   FUN = mean))

Fig2A <- ggplot(Poll_data, aes(x = Julian_day, y = Poll_spec, fill = Period)) + 
  geom_density_ridges(scale = 1, 
                      jittered_points = FALSE, 
                      point_fill = "grey80", 
                      alpha = .6,
                      quantile_lines = FALSE, quantiles = 2) +
  scale_fill_manual(values = c("grey70", "#D9af27")) +  
  labs(x = "Day of year", 
       y = "Pollinator species") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 

### (B) Point range plot ----

# Calculate mean, standard error, and 95% CI of Julian_day for each Species and Period
summary_data_poll <- Poll_data %>%
  group_by(Poll_spec, Period) %>%
  summarize(
    mean_julian_day = mean(Julian_day),
    sd_julian_day = sd(Julian_day),
    se_julian_day = sd_julian_day / sqrt(n()),
    ci_lower_mean = mean_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_mean = mean_julian_day + qt(0.975, df = n() - 1) * se_julian_day,
  )

# plot mean and 95% CI
Fig2B <- ggplot(summary_data_poll, aes(x = mean_julian_day,
                                       y = Poll_spec,
                                       color = Period,
                                       group = Period)) +
  geom_pointrange(aes(xmin = ci_lower_mean, xmax = ci_upper_mean),
                  position = position_dodge(0.5)) +
  scale_color_manual(values = c("grey70", "#D9af27")) +
  theme_minimal() +
  coord_cartesian(xlim = c(120, 260)) +
  labs(x = "Day of year") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

### (C) T-test results ----

## Statistical tests ----
t_test_results_poll <- Poll_data %>%
  group_by(Poll_spec) %>%
  summarize(
    p_value_poll = t.test(Julian_day ~ Period)$p.value,
    Stderr_poll = t.test(Julian_day ~ Period)$stderr,
    mean_past_poll = mean(Julian_day[Period == unique(Period)[1]]),
    mean_present_poll = mean(Julian_day[Period == unique(Period)[2]]),
    min_past_poll = min(Julian_day[Period == unique(Period)[1]]),
    min_present_poll = min(Julian_day[Period == unique(Period)[2]]),
    max_past_poll = max(Julian_day[Period == unique(Period)[1]]),
    max_present_poll = max(Julian_day[Period == unique(Period)[2]]),
    .groups = "drop") %>%
  mutate(
    significance = case_when(
      p_value_poll < 0.001 ~ "***",
      p_value_poll < 0.01 ~ "**",
      p_value_poll < 0.05 ~ "*",
      TRUE ~ "ns" # Not significant
    )
  )

# Visualize results
Fig2C <- ggplot(t_test_results_poll, aes(x = mean_present_poll - mean_past_poll, 
                                         y = Poll_spec)) +
  geom_point(aes(color = ifelse(p_value_poll < 0.05, "Significant", "Non-significant")), size = 2) +  
  geom_errorbar(aes(xmin = mean_present_poll - mean_past_poll - 1.96 * Stderr_poll,
                    xmax = mean_present_poll - mean_past_poll + 1.96 * Stderr_poll),
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
        axis.text.y = element_blank()) 

# Arrange plots 
Fig2A <- Fig2A + theme(legend.position = "none")
Fig2B <- Fig2B + theme(legend.position = "none")
Fig2C <- Fig2C + theme(legend.position = "none")

Fig2 <- cowplot::plot_grid(
  plot_grid(Fig2A, ncol = 1, rel_widths = c(5, 5)),  
  plot_grid(Fig2B, Fig2C, ncol = 2, rel_widths = c(0.1, 0.1)), 
  align = "hv",
  ncol = 2  ,
  label_size = 12,
  labels = "AUTO"
)

# ggsave("Fig_2.pdf", width = 8, height = 9, dpi = 300)

#_______________________________________________________________________

# 3: Shifts in plant phenology - individual species ----
#_______________________________________________________________________

## Fig. 3 ----

# Reorder Species by median Julian_day
Plant_data$plant_spec <- reorder(Plant_data$Plant_spec,
                                 ave(Plant_data$Julian_day,
                                     Plant_data$Plant_spec,
                                     FUN = mean))

### (A) Ridge plot ----
Fig3A <- ggplot(Plant_data, aes(x = Julian_day, y = plant_spec, fill = Period)) + 
  geom_density_ridges(scale = 1, 
                      jittered_points = FALSE, 
                      point_fill = "grey80", 
                      alpha = .6,
                      quantile_lines = FALSE, quantiles = 2) +
  scale_fill_manual(values = c("grey70", "springgreen4")) +  
  labs(x = "Day of year", 
       y = "Plant species") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") 


### (B) Point range plot ----

# Calculate mean, standard error, and 95% CI of Julian_day for each Species and Period
summary_data_plants <- Plant_data %>%
  group_by(Plant_spec, Period) %>%
  summarize(
    mean_julian_day = mean(Julian_day),
    sd_julian_day = sd(Julian_day),
    se_julian_day = sd_julian_day / sqrt(n()),
    ci_lower_mean = mean_julian_day - qt(0.975, df = n() - 1) * se_julian_day,
    ci_upper_mean = mean_julian_day + qt(0.975, df = n() - 1) * se_julian_day
  )

head(summary_data_plants)

custom_order <- c("Rubus deliciosus",
                  "Penstemon virens", "Penstemon secundiflorus",
                  "Rubus idaeus var. strigosus", "Penstemon glaber",
                  "Geranium richardsonii", "Aquilegia coerulea",
                  "Geranium caespitosum", "Aconitum columbianum",
                  "Chamaenerion angustifolium")


# Convert Species column to factor with custom order
summary_data_plants$Plant_spec <- factor(summary_data_plants$Plant_spec, levels = custom_order)

Fig3B <- ggplot(summary_data_plants, aes(x = mean_julian_day, 
                                         y = Plant_spec, 
                                         color = Period, 
                                         group = Period)) +
  geom_pointrange(aes(xmin = ci_lower_mean, xmax = ci_upper_mean),
                  position = position_dodge(0.5)) +
  scale_color_manual(values = c("grey70", "springgreen4")) + 
  theme_minimal() +
  coord_cartesian(xlim = c(120, 260)) +
  labs(x = "Day of year") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none", 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank()
  ) 

### (C) T-test results ----

## Statistical tests  ----
t_test_results_plants <- Plant_data %>%
  group_by(Plant_spec) %>%
  summarize(
    p_value_plants = t.test(Julian_day ~ Period)$p.value,
    Stderr_plants = t.test(Julian_day ~ Period)$stderr,
    mean_past_plants = mean(Julian_day[Period == unique(Period)[1]]),
    mean_present_plants = mean(Julian_day[Period == unique(Period)[2]]),
    min_past_plants = min(Julian_day[Period == unique(Period)[1]]),
    min_present_plants = min(Julian_day[Period == unique(Period)[2]]),
    max_past_plants = max(Julian_day[Period == unique(Period)[1]]),
    max_present_plants = max(Julian_day[Period == unique(Period)[2]]),
    .groups = "drop") %>%
  mutate(
    significance = case_when(
      p_value_plants < 0.001 ~ "***",
      p_value_plants < 0.01 ~ "**",
      p_value_plants < 0.05 ~ "*",
      TRUE ~ "ns" # Not significant
    )
  )

t_test_results_plants$Plant_spec <- factor(t_test_results_plants$Plant_spec, levels = custom_order)


# plot
Fig3C <- ggplot(t_test_results_plants, aes(x = mean_present_plants - mean_past_plants, 
                                           y = Plant_spec)) +
  geom_point(aes(color = ifelse(p_value_plants < 0.05, "Significant", "Non-significant")), size = 2) +
  geom_errorbar(aes(xmin = mean_present_plants - mean_past_plants - 1.96 * Stderr_plants,
                    xmax = mean_present_plants - mean_past_plants + 1.96 * Stderr_plants),
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
        axis.text.y = element_blank()
  ) 


# Arrange plots 

Fig3A <- Fig3A + theme(legend.position = "none")
Fig3B <- Fig3B + theme(legend.position = "none")
Fig3C <- Fig3C + theme(legend.position = "none")

Fig3 <- cowplot::plot_grid(
  plot_grid(Fig3A, ncol = 1, rel_widths = c(5, 5)),  
  plot_grid(Fig3B, Fig3C, ncol = 2, rel_widths = c(0.1, 0.1)), 
  align = "hv",
  ncol = 2  ,
  label_size = 12,
  labels = "AUTO"
)

# ggsave("Fig_3.pdf", width = 8, height = 9, dpi = 300)

#_______________________________________________________________________

# 4: Persistence of interactions ----
#_______________________________________________________________________

# make new variable "Interaction", concatenating Plant and Pollinator
Poll_data <- Poll_data %>% 
  unite("Interaction", Poll_spec, Plant_spec, sep = "-", remove = FALSE)

# calculate interaction persistence
int_persistence <- Poll_data %>%
  group_by(Interaction, Period) %>%
  summarise(
    Total_Freq = sum(Freq),
  ) %>%
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

# merge interaction persistence data and other phenological data
t_test_results_plants <- rename(t_test_results_plants, Plant_species = Plant_spec)
t_test_results_poll <- rename(t_test_results_poll, Poll_species = Poll_spec)

data <- merge(int_persistence, t_test_results_plants, by = "Plant_species")
data <- merge(data, t_test_results_poll, by = "Poll_species")
str(data)

# look at interactions that only occured in the past or the present
# Interactions not observed in the past but observed in the present
interaction_gains <- data %>%
  filter(freq_Past == 0 & freq_Present > 0)
dim(unique(interaction_gains)) # 28 interactions observed exclusively in the present

# Interactions observed in the past but not in the present
interaction_losses <- data %>%
  filter(freq_Past > 0 & freq_Present == 0)
dim(unique(interaction_losses)) # 49 interaction observed exclusively in the past

# calculate duration of activity period of plants and pollinators
data$poll_activity_period_past <- data$max_past_poll - data$min_past_poll
data$poll_activity_period_pres <- data$max_present_poll - data$min_present_poll
data$plant_activity_period_past <- data$max_past_plants - data$min_past_plants
data$plant_activity_period_pres <- data$max_present_plants - data$min_present_plants

# calculate days of overlap of plants and pollinators:
data$plant_poll_overlap_days_past <- pmax(0, pmin(data$max_past_plants, data$max_past_poll) - 
                                            pmax(data$min_past_plants, data$min_past_poll))
data$plant_poll_overlap_days_present <- pmax(0, pmin(data$max_present_plants, data$max_present_poll) - 
                                               pmax(data$min_present_plants, data$min_present_poll))
data$overlap_change <- data$plant_poll_overlap_days_present - data$plant_poll_overlap_days_past

# remove rows where interaction were not observed in the past
filtered_data <- data %>% filter(freq_Past > 0)

filtered_data$Plant_species <- as.factor(filtered_data$Plant_species)
filtered_data$Poll_species <- as.factor(filtered_data$Poll_species)


## Fit model ----
model <- glmer(interaction_persistence ~ overlap_change + (1|Plant_species) + (1|Poll_species), 
               data = filtered_data, 
               family = binomial)
summary(model)

## Fig. 4  ----
Fig4 <- ggplot(filtered_data, aes(x = overlap_change, y = interaction_persistence)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, color = "red") +
  labs(x = "Change in days of overlap", y = "interaction persistence") +
  xlim(-40, 70) +
  theme_minimal()

# ggsave("Fig_4.pdf", width = 5, height = 5, dpi = 300)

