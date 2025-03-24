
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
# This script contains code for reproducing the analyses and figures presented in the supplementary 
# material of the Manuscript:

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

#_______________________________________________________________________________________________

# Fig. S2 ----
## Density curves showing the distribution of observations for each day of the 
# year in the past and present.

FigS2 <- ggplot(Plant_data, aes(x = Julian_day, fill = Period)) +
  geom_density(alpha = 0.6, color = "white") +
  scale_fill_manual(values = c("Past" = "grey30", "Present" = "skyblue3")) +
  labs(
    x = "Day of year",
    y = "Density",
    fill = "Period"
  ) +
  xlim(c(140, 240)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

# ggsave("Fig_S2.pdf", width = 5, height = 6, dpi = 300)

# Fig. S3 ----
# Histogram showing the frequency of contemporary observations along 
# different elevations.

# read data
Pres <- read.csv("data/Pres.csv")

# histogram showing elevation of observation in present dataset
Pres$Ele_m <- as.numeric(Pres$Ele_m)
FigS3 <- hist(Pres$Ele_m, 
           breaks = c(1800, 1900, 2000, 2100, 2200, 2300, 2400,  
                      2500, 2600, 2700, 2800, 2900, 3000, 3100, 
                      3200, 3300, 3400, 3500, 3600, 3700, 3800, 
                      3900, 4000, 4100, 4200, 4300, 4400),
           xlab = "Elevation (m)",
           main = "")


# Fig. S4 ----
# Correlation of results obtained using the MinMax and percentile methods. 

## (A) Onset (Plants) ----
# Calculate correlation coefficients for periods seperately
(correlations <- Plant_metrics %>%
   group_by(Period) %>%
   summarize(
     Correlation = cor(Onset, Onset_percentile, method = "pearson")
   ) %>%
   mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4A <- ggplot(Plant_metrics, aes(x = Onset, y = Onset_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Onset (plants)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (B) Onset (Pollinators) ----
(correlations <- Poll_metrics %>%
   group_by(Period) %>%
   summarize(
     Correlation = cor(Onset, Onset_percentile, method = "pearson")
   ) %>%
   mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4B <- ggplot(Poll_metrics, aes(x = Onset, y = Onset_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Onset (pollinators)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (C) Central Tendency (Plants) ----
(correlations <- Plant_metrics %>%
    group_by(Period) %>%
    summarize(
      Correlation = cor(Mean, Median_percentile, method = "pearson")
    ) %>%
    mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4C <- ggplot(Plant_metrics, aes(x = Mean, y = Median_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Central tendency (plants)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (D) Central Tendency (Pollinators) ----
(correlations <- Poll_metrics %>%
   group_by(Period) %>%
   summarize(
     Correlation = cor(Mean, Median_percentile, method = "pearson")
   ) %>%
   mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4D <- ggplot(Poll_metrics, aes(x = Mean, y = Median_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Central tendency (pollinators)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (E) End (Plants) ----
(correlations <- Plant_metrics %>%
    group_by(Period) %>%
    summarize(
      Correlation = cor(End, End_percentile, method = "pearson")
    ) %>%
    mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4E <- ggplot(Plant_metrics, aes(x = End, y = End_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "End (plants)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (F) End (Pollinators) ----
(correlations <- Poll_metrics %>%
   group_by(Period) %>%
   summarize(
     Correlation = cor(End, End_percentile, method = "pearson")
   ) %>%
   mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4F <- ggplot(Poll_metrics, aes(x = End, y = End_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "End (pollinators)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(140, 240) +
  ylim(140, 240) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (G) Duration (Plants) ----
(correlations <- Plant_metrics %>%
    group_by(Period) %>%
    summarize(
      Correlation = cor(Duration, Duration_percentile, method = "pearson")
    ) %>%
    mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4G <- ggplot(Plant_metrics, aes(x = Duration, y = Duration_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Duration (plants)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(0, 105) +
  ylim(0, 90) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (H) Duration (Pollinators) ----
(correlations <- Poll_metrics %>%
    group_by(Period) %>%
    summarize(
      Correlation = cor(Duration, Duration_percentile, method = "pearson")
    ) %>%
    mutate(Label = paste0("r = ", round(Correlation, 2))))

FigS4H <- ggplot(Poll_metrics, aes(x = Duration, y = Duration_percentile, color = Period)) +
  geom_point(size = 3, aes(shape = Period), alpha = 1) + # Scatter points
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Period), size = 1) + # Correlation lines
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) + # Perfect correlation
  scale_color_manual(values = c("Past" = "grey50", "Present" = "skyblue3")) +
  scale_shape_manual(values = c("Past" = 16, "Present" = 17)) +
  scale_linetype_manual(values = c("Past" = "solid", "Present" = "solid")) +
  theme_minimal() +
  labs(title = "Duration (pollinators)",
       x = "MinMax method",
       y = "Percentile method",
       color = "Period",
       shape = "Period",
       linetype = "Period"
  ) +
  xlim(0, 105) +
  ylim(0, 90) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Arrange plots
FigS4A <- FigS4A + theme(legend.position = "none")
FigS4B <- FigS4B + theme(legend.position = "none")
FigS4C <- FigS4C + theme(legend.position = "none")
FigS4D <- FigS4D + theme(legend.position = "none")
FigS4E <- FigS4E + theme(legend.position = "none")
FigS4F <- FigS4F + theme(legend.position = "none")
FigS4G <- FigS4G + theme(legend.position = "none")
FigS4H <- FigS4H + theme(legend.position = "none")

Fig_S4 <- cowplot::plot_grid(FigS4B, FigS4A, FigS4D, FigS4C,
                             FigS4F, FigS4E, FigS4H, FigS4G,
                             labels = "AUTO",
                             ncol = 2, nrow = 4,
                             label_size = 12)

# ggsave("Fig_S4.pdf", width = 6, height = 9, dpi = 300)


# Fig. S5 ----
# Point plots showing the means and standard deviations of onset, mean, end and 
# duration of pollinator flight and plant flowering, using the percentile approach

## (A) Onset ----
FigS5A <- ggplot(metrics_mean_se, aes(x = Period, y = mean_onset_perc, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_onset_perc - se_onset_perc, ymax = mean_onset_perc + se_onset_perc), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "Onset (5th percentile)") +
  theme_minimal() +
  ylim(140, 250) 

## (B) Median ----
FigS5B <- ggplot(metrics_mean_se, aes(x = Period, y = mean_median_perc, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_median_perc - se_median_perc, ymax = mean_median_perc + se_median_perc), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "Median (50th percentile)") +
  theme_minimal() +
  ylim(140, 250) 

## (C) End ----
FigS5C <- ggplot(metrics_mean_se, aes(x = Period, y = mean_end_perc, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_end_perc - se_end_perc, ymax = mean_end_perc + se_end_perc), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "End (95th percentile)") +
  theme_minimal()+
  ylim(140, 250) 

## (D) Duration ----
FigS5D <- ggplot(metrics_mean_se, aes(x = Period, y = mean_Duration_percentile, color = Type, shape = Type, group = Type)) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_Duration_percentile - se_Duration_percentile, ymax = mean_Duration_percentile + se_Duration_percentile), 
                width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = c("Plants" = "springgreen4", "Pollinators" = "#D9af27")) +
  scale_shape_manual(values = c("Plants" = 17, "Pollinators" = 16)) + 
  labs(x = "",
       y = "95th - 5th percentile") +
  theme_minimal() +
  theme(legend.position = "right") +
  ylim(19, 65)

# Arrange plots
FigS5A <- FigS5A + theme(legend.position = "none")
FigS5B <- FigS5B + theme(legend.position = "none")
FigS5C <- FigS5C + theme(legend.position = "none")
FigS5D <- FigS5D + theme(legend.position = "none")

FigS5 <- cowplot::plot_grid(FigS5A, FigS5B , FigS5C, FigS5D, 
                            labels = "AUTO",
                            ncol = 2, nrow = 2,
                            label_size = 12)

#ggsave("Fig_S5.pdf", width = 5, height = 6, dpi = 300)

## Statistical tests ----
past_data <- Metrics %>% filter(Period == "Past")
present_data <- Metrics %>% filter(Period == "Present")
poll_data <- Metrics %>% filter(Type == "Pollinators")
plant_data <- Metrics %>% filter(Type == "Plants")

{
  percentile_test_results <- list()
  # onset
  (percentile_test_results$t_test_onset_past <- t.test(Onset_percentile ~ Type, data = past_data))
  (percentile_test_results$t_test_onset_poll <- t.test(Onset_percentile ~ Period, data = poll_data))
  (percentile_test_results$t_test_onset_present <- t.test(Onset_percentile ~ Type, data = present_data))
  (percentile_test_results$t_test_onset_plant <- t.test(Onset_percentile ~ Period, data = plant_data))
  # median
  (percentile_test_results$t_test_median_past <- t.test(Median_percentile ~ Type, data = past_data))
  (percentile_test_results$t_test_median_poll <- t.test(Median_percentile ~ Period, data = poll_data))
  (percentile_test_results$t_test_median_present <- t.test(Median_percentile ~ Type, data = present_data))
  (percentile_test_results$t_test_median_plant <- t.test(Median_percentile ~ Period, data = plant_data))
  # end
  (percentile_test_results$t_test_end_past <- t.test(End_percentile ~ Type, data = past_data))
  (percentile_test_results$t_test_end_poll <- t.test(End_percentile ~ Period, data = poll_data))
  (percentile_test_results$t_test_end_present <- t.test(End_percentile ~ Type, data = present_data))
  (percentile_test_results$t_test_end_plant <- t.test(End_percentile ~ Period, data = plant_data))
  # duration
  (percentile_test_results$t_test_duration_past <- t.test(Duration_percentile ~ Type, data = past_data))
  (percentile_test_results$t_test_duration_poll <- t.test(Duration_percentile ~ Period, data = poll_data))
  (percentile_test_results$t_test_duration_present <- t.test(Duration_percentile ~ Type, data = present_data))
  (percentile_test_results$t_test_duration_plant <- t.test(Duration_percentile ~ Period, data = plant_data))
  
  percentile_test_results
  }

# Fig. S6 ----
# Correlations between phenological metrics of plants and pollinators in the 
# past and in the present using both the MinMax and the percentile approach. 

## (A) Onset (MinMax) ----
# Reshape data to have Past and Present onset as separate columns
  Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, Onset) %>%
    pivot_wider(names_from = Period, values_from = Onset) %>%
    rename(Past_Onset = Past, Present_Onset = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, Onset) %>%
    pivot_wider(names_from = Period, values_from = Onset) %>%
    rename(Past_Onset = Past, Present_Onset = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)

FigS6A <- ggplot(Combined_metrics, aes(x = Past_Onset, 
                                       y = Present_Onset, 
                                       color = Type, 
                                       shape = Type)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
  scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
  scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
  theme_minimal() +
  labs(title = "MinMax method",
       x = "Onset past",
       y = "Onset present",
       color = "Type",
       shape = "Type",
       linetype = "Type"
  ) +
  xlim(min(Combined_metrics$Past_Onset, na.rm = TRUE) - 5, 
       max(Combined_metrics$Past_Onset, na.rm = TRUE) + 5) +
  ylim(min(Combined_metrics$Present_Onset, na.rm = TRUE) - 5, 
       max(Combined_metrics$Present_Onset, na.rm = TRUE) + 5) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (B) Mean (MinMax) ----
  Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, Mean) %>%
    pivot_wider(names_from = Period, values_from = Mean) %>%
    rename(Past_Mean = Past, Present_Mean = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, Mean) %>%
    pivot_wider(names_from = Period, values_from = Mean) %>%
    rename(Past_Mean = Past, Present_Mean = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)

FigS6B <- ggplot(Combined_metrics, aes(x = Past_Mean, 
                                       y = Present_Mean, 
                                       color = Type, 
                                       shape = Type)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
  scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
  scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
  theme_minimal() +
  labs(x = "Central tendency past",
       y = "Central tendency present",
       color = "Type",
       shape = "Type",
       linetype = "Type"
  ) +
  xlim(min(Combined_metrics$Past_Mean, na.rm = TRUE) - 5, 
       max(Combined_metrics$Past_Mean, na.rm = TRUE) + 5) +
  ylim(min(Combined_metrics$Present_Mean, na.rm = TRUE) - 5, 
       max(Combined_metrics$Present_Mean, na.rm = TRUE) + 5) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (C) End (MinMax) ----
Poll_metrics_wide <- Poll_metrics %>%
  dplyr::select(Species, Period, End) %>%
  pivot_wider(names_from = Period, values_from = End) %>%
  rename(Past_End = Past, Present_End = Present)

Plant_metrics_wide <- Plant_metrics %>%
  dplyr::select(Species, Period, End) %>%
  pivot_wider(names_from = Period, values_from = End) %>%
  rename(Past_End = Past, Present_End = Present)

Poll_metrics_wide$Type <- "Pollinator"
Plant_metrics_wide$Type <- "Plant"

Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)

FigS6C <- ggplot(Combined_metrics, aes(x = Past_End, y = Present_End, color = Type, shape = Type)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
  scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
  scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
  theme_minimal() +
  labs(x = "End past",
       y = "End present",
       color = "Type",
       shape = "Type",
       linetype = "Type"
  ) +
  xlim(min(Combined_metrics$Past_End, na.rm = TRUE) - 5, 
       max(Combined_metrics$Past_End, na.rm = TRUE) + 5) +
  ylim(min(Combined_metrics$Present_End, na.rm = TRUE) - 5, 
       max(Combined_metrics$Present_End, na.rm = TRUE) + 5) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (D) Duration (MinMax) ----
  Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, Duration) %>%
    pivot_wider(names_from = Period, values_from = Duration) %>%
    rename(Past_Duration = Past, Present_Duration = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, Duration) %>%
    pivot_wider(names_from = Period, values_from = Duration) %>%
    rename(Past_Duration = Past, Present_Duration = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)

FigS6D <- ggplot(Combined_metrics, aes(x = Past_Duration, 
                                       y = Present_Duration, 
                                       color = Type, 
                                       shape = Type)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
  scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
  scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
  theme_minimal() +
  labs(x = "Duration past",
       y = "Duration present",
       color = "Type",
       shape = "Type",
       linetype = "Type"
  ) +
  xlim(min(Combined_metrics$Past_Duration, na.rm = TRUE) - 5, 
       max(Combined_metrics$Past_Duration, na.rm = TRUE) + 5) +
  ylim(min(Combined_metrics$Present_Duration, na.rm = TRUE) - 5, 
       max(Combined_metrics$Present_Duration, na.rm = TRUE) + 5) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

## (E) Onset (Percentile) ----
Poll_metrics_wide <- Poll_metrics %>%
  dplyr::select(Species, Period, Onset_percentile) %>%
  pivot_wider(names_from = Period, values_from = Onset_percentile) %>%
  rename(Past_Onset_perc = Past, Present_Onset_perc = Present)

Plant_metrics_wide <- Plant_metrics %>%
  dplyr::select(Species, Period, Onset_percentile) %>%
  pivot_wider(names_from = Period, values_from = Onset_percentile) %>%
  rename(Past_Onset_perc = Past, Present_Onset_perc = Present)

Poll_metrics_wide$Type <- "Pollinator"
Plant_metrics_wide$Type <- "Plant"

Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)

FigS6E <- ggplot(Combined_metrics, aes(x = Past_Onset_perc, 
                                       y = Present_Onset_perc, 
                                       color = Type, 
                                       shape = Type)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
  scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
  scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
  theme_minimal() +
  labs(title = "Percentile method",
       x = "Onset past",
       y = "Onset present",
       color = "Type",
       shape = "Type",
       linetype = "Type"
  ) +
  xlim(min(Combined_metrics$Past_Onset_perc, na.rm = TRUE) - 5, 
       max(Combined_metrics$Past_Onset_perc, na.rm = TRUE) + 5) +
  ylim(min(Combined_metrics$Present_Onset_perc, na.rm = TRUE) - 5, 
       max(Combined_metrics$Present_Onset_perc, na.rm = TRUE) + 5) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )


## (F) Median (Percentile) ----
Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, Median_percentile) %>%
    pivot_wider(names_from = Period, values_from = Median_percentile) %>%
    rename(Past_Median = Past, Present_Median = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, Median_percentile) %>%
    pivot_wider(names_from = Period, values_from = Median_percentile) %>%
    rename(Past_Median = Past, Present_Median = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)
  
  # plot
  FigS6F <- ggplot(Combined_metrics, aes(x = Past_Median, 
                                         y = Present_Median, 
                                         color = Type, 
                                         shape = Type)) +
    geom_point(size = 3, alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
    scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
    scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
    scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
    theme_minimal() +
    labs(x = "Central tendency past",
         y = "Central tendency present",
         color = "Type",
         shape = "Type",
         linetype = "Type"
    ) +
    xlim(min(Combined_metrics$Past_Median, na.rm = TRUE) - 5, 
         max(Combined_metrics$Past_Median, na.rm = TRUE) + 5) +
    ylim(min(Combined_metrics$Present_Median, na.rm = TRUE) - 5, 
         max(Combined_metrics$Present_Median, na.rm = TRUE) + 5) +
    theme(
      legend.position = "top",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )

## (G) End (Percentile) ----
Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, End_percentile) %>%
    pivot_wider(names_from = Period, values_from = End_percentile) %>%
    rename(Past_End_perc = Past, Present_End_perc = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, End_percentile) %>%
    pivot_wider(names_from = Period, values_from = End_percentile) %>%
    rename(Past_End_perc = Past, Present_End_perc = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)
  
  FigS6G <- ggplot(Combined_metrics, aes(x = Past_End_perc, 
                                         y = Present_End_perc, 
                                         color = Type, 
                                         shape = Type)) +
    geom_point(size = 3, alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
    scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
    scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
    scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
    theme_minimal() +
    labs(x = "End past",
         y = "End present",
         color = "Type",
         shape = "Type",
         linetype = "Type"
    ) +
    xlim(min(Combined_metrics$Past_End_perc, na.rm = TRUE) - 5, 
         max(Combined_metrics$Past_End_perc, na.rm = TRUE) + 5) +
    ylim(min(Combined_metrics$Present_End_perc, na.rm = TRUE) - 5, 
         max(Combined_metrics$Present_End_perc, na.rm = TRUE) + 5) +
    theme(
      legend.position = "top",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )

## (H) Duration (Percentile) ----
Poll_metrics_wide <- Poll_metrics %>%
    dplyr::select(Species, Period, Duration_percentile) %>%
    pivot_wider(names_from = Period, values_from = Duration_percentile) %>%
    rename(Past_Duration_perc = Past, Present_Duration_perc = Present)
  
  Plant_metrics_wide <- Plant_metrics %>%
    dplyr::select(Species, Period, Duration_percentile) %>%
    pivot_wider(names_from = Period, values_from = Duration_percentile) %>%
    rename(Past_Duration_perc = Past, Present_Duration_perc = Present)
  
  Poll_metrics_wide$Type <- "Pollinator"
  Plant_metrics_wide$Type <- "Plant"
  
  # Combine the datasets
  Combined_metrics <- rbind(Poll_metrics_wide, Plant_metrics_wide)
  
  # Make the combined scatter plot
  FigS6H <- ggplot(Combined_metrics, aes(x = Past_Duration_perc, 
                                         y = Present_Duration_perc, 
                                         color = Type, 
                                         shape = Type)) +
    geom_point(size = 3, alpha = 1) +
    geom_smooth(method = "lm", se = FALSE, aes(linetype = Type), size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.8) +
    scale_color_manual(values = c("Pollinator" = "#D9af27", "Plant" = "springgreen4")) +
    scale_shape_manual(values = c("Pollinator" = 16, "Plant" = 17)) +
    scale_linetype_manual(values = c("Pollinator" = "solid", "Plant" = "solid")) +
    theme_minimal() +
    labs(x = "Duration past",
         y = "Duration present",
         color = "Type",
         shape = "Type",
         linetype = "Type"
    ) +
    xlim(min(Combined_metrics$Past_Duration_perc, na.rm = TRUE) - 5, 
         max(Combined_metrics$Past_Duration_perc, na.rm = TRUE) + 5) +
    ylim(min(Combined_metrics$Present_Duration_perc, na.rm = TRUE) - 5, 
         max(Combined_metrics$Present_Duration_perc, na.rm = TRUE) + 5) +
    theme(
      legend.position = "top",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )

# Arrange plots
FigS6A <- FigS6A + theme(legend.position = "none")
FigS6B <- FigS6B + theme(legend.position = "none")
FigS6C <- FigS6C + theme(legend.position = "none")
FigS6D <- FigS6D + theme(legend.position = "none")
FigS6E <- FigS6E + theme(legend.position = "none")
FigS6F <- FigS6F + theme(legend.position = "none")
FigS6G <- FigS6G + theme(legend.position = "none")
FigS6H <- FigS6H + theme(legend.position = "none")
  
Fig_S6 <- cowplot::plot_grid(FigS6A, FigS6E, FigS6B, FigS6F,
                             FigS6C, FigS6G, FigS6D, FigS6H,
                             labels = "AUTO",
                             ncol = 2, nrow = 4,
                             label_size = 12)

# ggsave("Fig_S6.pdf", width = 8, height = 11, dpi = 300)

