
############.
## Set up ## -----
############.

rm(list = ls())

setwd("/home/ubuntu/R_code")
# setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Code/Scenario-hub-visualisations")

# library(ggplot2)
# library(dplyr)
# library(tidyr)
library(tidyverse)
library(patchwork)
library(lubridate)
library(reshape2)
library(tictoc)

# Read in data from another script
source('Process_data.R') # Read in simulation results from both teams
source("Colour_palette.R")  # Define colour palettes and visualisation standards for outputs

# ################################################################################.
# # Adding reductions in infections, cases and hospitalisations to the dataframe  # ------ 
# ################################################################################.

add_reductions <- function(df){
  df <- df |>
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
  
  return(df)
}

dat_uom_list <- lapply(dat_uom_list, add_reductions)
dat_uom_inf <- dat_uom_list[[1]] %>%
  mutate(team = "no waning")
dat_uom_dis <- dat_uom_list[[2]] %>%
  mutate(team = "no waning")
dat_uom_adm <- dat_uom_list[[3]] %>%
  mutate(team = "no waning")
dat_uom_fatal <- dat_uom_list[[4]] %>%
  mutate(team = "no waning")

dat_uom_list_with_waning <- lapply(dat_uom_list_with_waning, add_reductions)
dat_uom_inf_with_waning <- dat_uom_list_with_waning[[1]] %>%
  mutate(team = "with waning")
dat_uom_dis_with_waning <- dat_uom_list_with_waning[[2]] %>%
  mutate(team = "with waning")
dat_uom_adm_with_waning <- dat_uom_list_with_waning[[3]] %>%
  mutate(team = "with waning")
dat_uom_fatal_with_waning <- dat_uom_list_with_waning[[4]] %>%
  mutate(team = "with waning")

dat_uom_inf <- dat_uom_inf %>%
  bind_rows(dat_uom_inf_with_waning)

dat_uom_dis <- dat_uom_dis %>%
  bind_rows(dat_uom_dis_with_waning)

dat_uom_adm <- dat_uom_adm %>%
  bind_rows(dat_uom_adm_with_waning)

dat_uom_fatal <- dat_uom_fatal %>%
  bind_rows(dat_uom_fatal_with_waning)

dat_uom_list <- list(dat_uom_inf, dat_uom_dis, dat_uom_adm, dat_uom_fatal)

###########.
## Plots  ## -----
###########.

dir.create("Plots", showWarnings = FALSE)

y_axis_name <- c(
  "Infections",
  "Cases",
  "Hospitalisations",
  "Fatalities"
)

for ( i in seq_along(dat_uom_list)){
  
  fig_data <- dat_uom_list[[i]] |>
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
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = y_axis_name[i], title = paste(y_axis_name[i], "across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = FALSE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste(y_axis_name[i], ".png", sep = "")), height = 10, width = 14, dpi = 300)
  
  fig_data <- fig_data |>
    filter(scenario != "Status_quo")
  
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
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = paste(y_axis_name[i], "averted"), title = paste(y_axis_name[i], "averted across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = TRUE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste(y_axis_name[i], "_averted.png", sep = "")), height = 10, width = 14, dpi = 300)
  
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
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = "Age groups", y = paste("%", y_axis_name[i], "averted"), title = paste( "%",y_axis_name[i], "averted across all age groups")) +
    facet_wrap(~scenario, ncol = 3, drop = TRUE, labeller = labeller(scenario = labeller_scenario))
  
  ggsave(file.path("Plots", paste("perc_",y_axis_name[i], "_averted.png", sep = "")), height = 10, width = 14, dpi = 300)
  
}

cat("Plots production finished.\n\n")
