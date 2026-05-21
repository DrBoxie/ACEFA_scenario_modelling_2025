
############.
## Set up ## -----
############.

rm(list = ls())

setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Scenario-hub-visualisations")

library(ggplot2)
#library(tidyverse)
library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)
library(socialmixr)
library(reshape2)
library(tidyverse)
library(tictoc)

# Read in data from another script
source('Process data.R') # Read in simulation results from both teams
source("Colour_palette.R")  # Define colour palettes and visualisation standards for outputs

# ################################################################################.
# # Adding reductions in infections, cases and hospitalisations to the dataframe  # ------ 
# ################################################################################.

dat_uom <- dat_uom |>
  group_by(simulation_index, age_group, target, team) |>
  mutate(
    reduction = if_else(
      scenario == "Status_quo",
      NA_real_,
      value[scenario == "Status_quo"] - value
    )
  ,
  perc_reduction = if_else(
    scenario == "Status_quo",
    NA_real_,
    100 * (value[scenario == "Status_quo"] - value) / value[scenario == "Status_quo"]
  )
  ) |>
  ungroup()

###########.
## Plots  ## -----
###########.

for (t in c("infection_incidence", "disease_incidence", "admission_incidence")){

# for (t in c("infection_incidence")){
  y_axis_name <- switch(t,
                        "infection_incidence" = "Infections",
                        "disease_incidence" = "Cases",
                        "admission_incidence" = "Hospitalisations"
                        )
  
  fig_data <- dat_uom |>
    filter(target == t) |>
    mutate(scenario = factor(scenario, levels = order_scenario2))
  
  # Absolute numbers for each target
  ggplot(fig_data) +
    geom_point(aes(x = age_group, y = value, colour = team), position = position_jitter(width = 0.1, height = 0), alpha = 0.02) +
    stat_summary(
      aes(x = age_group, y = value, colour = team),
      fun = median,
      geom = "point",
      shape = 95,
      size = 10
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3, alpha = 0.7) +
    scale_colour_manual(values = scenario_palettes$team) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect("white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = y_axis_name, title = paste(y_axis_name, "across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = FALSE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste(y_axis_name, ".png", sep = "")), height = 10, width = 14, dpi = 300)
  
  fig_data <- fig_data |>
    filter(!scenario == "Status_quo")
  
  # Reduction for each target
  ggplot(fig_data) +
    geom_point(aes(x = age_group, y = reduction, colour = team), position = position_jitter(width = 0.1, height = 0), alpha = 0.02) +
    stat_summary(
      aes(x = age_group, y = reduction, colour = team),
      fun = median,
      geom = "point",
      shape = 95,
      size = 10
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3, alpha = 0.7) +
    scale_colour_manual(values = scenario_palettes$team) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect("white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = paste(y_axis_name, "averted"), title = paste(y_axis_name, "averted across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = TRUE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste(y_axis_name, "_averted.png", sep = "")), height = 10, width = 14, dpi = 300)
  
  # Percentage reduction for each target
  ggplot(fig_data) +
    geom_point(aes(x = age_group, y = perc_reduction, colour = team), position = position_jitter(width = 0.1, height = 0), alpha = 0.02) +
    stat_summary(
      aes(x = age_group, y = perc_reduction, colour = team),
      fun = median,
      geom = "point",
      shape = 95,
      size = 10
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3, alpha = 0.7) +
    scale_colour_manual(values = scenario_palettes$team) +
    scale_y_continuous(labels = scales::comma) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect("white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = paste("%", y_axis_name, "averted"), title = paste( "%",y_axis_name, "averted across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = TRUE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste("perc_",y_axis_name, "_averted.png", sep = "")), height = 10, width = 14, dpi = 300)
  
}
