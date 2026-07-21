rm(list = ls())

source_file_type <- "with waning"

target_list <- c("infection", "disease", "admission")

library(ggplot2)
library(patchwork)
library(arrow)
library(RColorBrewer)
library(dplyr)

source_file_formatted <- gsub("\\s+", "_", source_file_type)

for (target in target_list){
  
  if (.Platform$OS.type == "windows"){
    setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Code")
    output_root <- file.path(getwd(), "Scenario-hub-visualisations", "data")
    source_file_path <- file.path(paste("Forward projection",source_file_type), paste("dat_uom_", source_file_formatted, sep = "_"))
  } else{
    setwd("/home/ubuntu/R_code")
    ROOT <- "/pvol"
    source_file_path <- file.path(ROOT, paste("R outputs", source_file_type), paste("dat_uom",target, source_file_formatted, sep="_"))
    output_path <- file.path(ROOT, "R code Plots")
  }
  
  # Read in the modelled outputs
  
  df <- open_dataset(source_file_path) %>%
    filter(scenario == "Status_quo") %>%
    select(-c(setting, team, scenario, target)) %>%
    group_by(simulation_index, run_nr, horizon) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    collect()
  
  df <- df %>%
    arrange(simulation_index, run_nr, horizon) %>%
    group_by(simulation_index, run_nr) %>%
    mutate(cum_val = cumsum(value)) %>%
    ungroup()
  
  if (target == "infection") {
    df_inf <- df
  } else if (target == "disease") {
    df_dis <- df
  } else if (target == "admission") {
    df_adm <- df
  }
}

rm(df)

p_inf <- ggplot(
  df_inf[df_inf$run_nr == 0, ],
  aes(x = horizon, y = value, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  coord_cartesian(ylim = c(0, 14000)) +
  ggtitle("Infections")

p_dis <- ggplot(
  df_dis[df_dis$run_nr == 0, ],
  aes(x = horizon, y = value, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  coord_cartesian(ylim = c(0, 400)) +
  ggtitle("Cases")

p_adm <- ggplot(
  df_adm[df_adm$run_nr == 0, ],
  aes(x = horizon, y = value, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle("Hospitalisations")

p_inf + p_dis + p_adm

ggsave(file.path(output_path, paste("baseline_infs_cases_hosps.png", sep = "")), height = 10, width = 14, dpi = 300)

p_cumul_inf <- ggplot(
  df_inf[df_inf$run_nr == 0, ],
  aes(x = horizon, y = cum_val, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  ggtitle("Infections")

p_cumul_dis <- ggplot(
  df_dis[df_dis$run_nr == 0, ],
  aes(x = horizon, y = cum_val, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  ggtitle("Cases")

p_cumul_adm <- ggplot(
  df_adm[df_adm$run_nr == 0, ],
  aes(x = horizon, y = cum_val, group = simulation_index)
) +
  geom_line(alpha = 0.2, colour = "black") +
  ggtitle("Hospitalisations")

p_cumul_inf + p_cumul_dis + p_cumul_adm


ggsave(file.path(output_path, paste("baseline_cumulative_infs_cases_hosps.png", sep = "")), height = 10, width = 14, dpi = 300)

##############################################################################################################
