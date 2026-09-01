rm(list = ls())

source_file_type <- "with waning"
vaccination_term <- "term1 vaccination"

target_list <- c("infection", "admission", "admission_best_flucan", "admission_worst_flucan")

library(ggplot2)
library(patchwork)
library(arrow)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

source_file_formatted <- gsub("\\s+", "_", source_file_type)

age_levels <- c("<1", "1-4", "5-11", "12-17", "18-64", "65-79", "80+", "All")

for (target in target_list){
  
  if (.Platform$OS.type == "windows"){
    setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Code")
    output_root <- file.path(getwd(), "Scenario-hub-visualisations", "data")
    source_file_path <- file.path(paste("Forward projection",source_file_type), paste("df_incid_", source_file_formatted, ".parquet", sep = ""))
  } else{
    setwd("/home/ubuntu/R_code")
    ROOT <- "/pvol"
    source_file_path <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), paste(target, "ar_ds", sep="_"))
    output_path <- file.path(ROOT, "R code Plots")
  }
  
  # Read in the modelled outputs
  
  df <- open_dataset(source_file_path) %>%
    collect()
  
  df_all_ages <- df %>%
    group_by(scenario, simulation_index, run_nr) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(age_group = "All")
  
  df <- dplyr::bind_rows(df, df_all_ages)
  rm(df_all_ages)
  
  df$age_group <- factor(df$age_group, levels = age_levels)
  
  if (target == "infection") {
    df_inf <- df
  } else if (target == "disease") {
    df_dis <- df
  } else if (target == "admission") {
    df_adm <- df
  } else if (target == "admission_best_flucan") {
    df_adm_best_flucan <- df
  } else if (target == "admission_worst_flucan") {
    df_adm_worst_flucan <- df
  }
}

rm(df)
invisible(gc())

particles_file_path <- file.path(ROOT, paste("Forward projection with waning", vaccination_term), "particles.csv")
particles <- read.csv(particles_file_path) %>%
  select(-X) %>%
  rownames_to_column(var = "simulation_index") %>%
  mutate(simulation_index = as.integer(simulation_index)-1)

ages <- read_csv(file.path("data", "age_pops.csv")) %>%
  mutate(`1-4` = `1` + `2-4`) %>%
  select(-`1`, -`2-4`) %>%
  pivot_longer(cols = everything(), names_to = "age_group",
               values_to = "pop_size") %>%
  bind_rows(
    summarise(
      .,
      age_group = "All",
      across(where(is.numeric), sum)
    )
  )

##############################################################################################################

##############################################################################################################

# Consistent text sizes across all figures
text_size_axis_title <- 17
text_size_axis_text  <- 15
text_size_other      <- 15
text_size_panel_tag  <- 20

##############################################################################################################

#################################################################################################################
  
## Manually creating lookup table for setting labels for different scenarios
# These are the actual labels being used
scenarios_considered <- c("Status_quo",
                          "Pessimistic_cover20_LAIV_5-12yo",
                          "Pessimistic_cover40_LAIV_5-12yo",
                          "Pessimistic_cover60_LAIV_5-12yo",
                          "Pessimistic_cover80_LAIV_5-12yo",
                          "Pessimistic_cover20_LAIV_5-18yo",
                          "Pessimistic_cover40_LAIV_5-18yo",
                          "Pessimistic_cover60_LAIV_5-18yo",
                          "Pessimistic_cover80_LAIV_5-18yo",
                          "Central_cover20_LAIV_5-12yo",
                          "Central_cover40_LAIV_5-12yo",
                          "Central_cover60_LAIV_5-12yo",
                          "Central_cover80_LAIV_5-12yo",
                          "Central_cover20_LAIV_5-18yo",
                          "Central_cover40_LAIV_5-18yo",
                          "Central_cover60_LAIV_5-18yo",
                          "Central_cover80_LAIV_5-18yo",
                          "Optimistic_cover20_LAIV_5-12yo",
                          "Optimistic_cover40_LAIV_5-12yo",
                          "Optimistic_cover60_LAIV_5-12yo",
                          "Optimistic_cover80_LAIV_5-12yo",
                          "Optimistic_cover20_LAIV_5-18yo",
                          "Optimistic_cover40_LAIV_5-18yo",
                          "Optimistic_cover60_LAIV_5-18yo",
                          "Optimistic_cover80_LAIV_5-18yo",
                          "Central_cover40_LAIV_2-5yo"
                          )

# These will be the new labels used
scen_age <- c("Baseline",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-12 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 5-18 year olds",
              "LAIV in 2-5 year olds"
              )

scen_effect <- c("Baseline",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Pessimistic",
                 "Central",
                 "Central",
                 "Central",
                 "Central",
                 "Central",
                 "Central",
                 "Central",
                 "Central",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Optimistic",
                 "Central"
                 )

scen_coverage <- c("Baseline",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "20%",
                   "40%",
                   "60%",
                   "80%",
                   "40%"
                   )

lookup_table <- data.frame(scenario = scenarios_considered,
                           scen_age = scen_age,
                           scen_effect = scen_effect,
                           scen_coverage = scen_coverage)

lookup_table$scen_effect <- factor(lookup_table$scen_effect,
                                   levels=c("Baseline",
                                            "Pessimistic",
                                            "Central",
                                            "Optimistic"))

lookup_table_half_life <- particles %>%
  select(half_life) %>%
  mutate(simulation_index = row_number()-1)%>%
  mutate(
    half_life_category = case_when(
      half_life <= 730  ~ "<= 2y",
      half_life > 730 & half_life <= 1095 ~ "2y - 3y",
      half_life > 1095 & half_life <= 1460 ~ "3y - 4y",
      half_life > 1460 ~ "> 4y"
    )
  ) %>%
  select(-half_life)

lookup_table_AR <- df_inf %>%
  filter(age_group == "All") %>%
  group_by(scenario, simulation_index) %>%
  summarise(AR = median(value), .groups = "drop")  %>%
  filter(scenario == "Status_quo")%>%
  mutate(AR_category = case_when(
    AR >= 100000 & AR < 250000 ~ "10-25",
    AR >= 250000 & AR < 350000 ~ "25-35",
    AR >= 350000               ~ "35-45" 
  )) %>%
  select(simulation_index, AR_category)

# Strength of seasonality category is each 1/3 of the range of beta_1
# Same for beta_0_x_beta_1
beta1_min <- min(particles$beta_1)
beta1_max <- max(particles$beta_1)
beta1_bin_size <- (beta1_max - beta1_min)/3
weak_beta1_max <- beta1_min + beta1_bin_size
medium_beta1_max <- weak_beta1_max + beta1_bin_size

beta0_x_beta1_min <- beta1_min * min(particles$beta_0)
beta0_x_beta1_max <- beta1_max * max(particles$beta_0)
beta0_x_beta1_bin_size <- (beta0_x_beta1_max - beta0_x_beta1_min)/3
weak_beta0_x_beta1_max <- beta0_x_beta1_min + beta0_x_beta1_bin_size
medium_beta0_x_beta1_max <- weak_beta0_x_beta1_max + beta0_x_beta1_bin_size

lookup_table_seasonality <- particles %>%
  select(simulation_index, beta_0, beta_1) %>%
  mutate(
    beta_0_x_beta_1 = beta_0 * beta_1,
    beta_1_category = case_when(
      beta_1 < weak_beta1_max ~ "weak",
      beta_1 >= weak_beta1_max & beta_1 < medium_beta1_max ~ "medium",
      beta_1 >= medium_beta1_max ~ "strong"
    ),
    beta0_x_beta1_category = case_when(
      beta_0_x_beta_1 < weak_beta0_x_beta1_max ~ "weak",
      beta_0_x_beta_1 >= weak_beta0_x_beta1_max & beta_0_x_beta_1 < medium_beta0_x_beta1_max ~ "medium",
      beta_0_x_beta_1 >= medium_beta0_x_beta1_max ~ "strong"
    )
    ) %>%
  select(-c(beta_0, beta_1, beta_0_x_beta_1))

t_vals <- seq(0, 364)
beta_peak_timing <- particles %>%
  crossing(t = t_vals)%>%
  mutate(
    beta = sin(2 * pi * (t - shift) / 365)
  ) %>%
  select(-c(beta_0, beta_1, shift, frac_S1, half_life, phi_S1, phi_V0, phi_V1)) %>%
  group_by(simulation_index) %>%
  slice_max(
    order_by = beta,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(-beta)

lookup_table_peak_beta <- beta_peak_timing %>%
  mutate(peak_beta_category = case_when(
    t < 180 ~ "< 180",
    t >= 180 & t < 210 ~ "180-210",
    t >= 210 & t < 240 ~ "210-240",
    t >= 240 ~ ">= 240",
  )
  ) %>%
  select(simulation_index, peak_beta_category)

################################################################################

################################################################################
## Summarising

for (target in target_list) {

  cat("Loading", target, "data.\n")
    
  if (target == "infection") {
    df <- df_inf
  } else if (target == "disease") {
    df <- df_dis
  } else if (target == "admission") {
    df <- df_adm
  } else if (target == "admission_best_flucan") {
    df <- df_adm_best_flucan
  } else if (target == "admission_worst_flucan") {
    df <- df_adm_worst_flucan
  }
  
  df <- df %>%
    left_join(ages, by = "age_group") %>%
    mutate(value = 100000 * value / pop_size) %>%
    select(-pop_size)
  
  # Seperate dataframes for baseline vs scenarios
  df_B <- df[df$scenario=="Status_quo",]
  df_S <- df[df$scenario!="Status_quo",]  
  
  # Join data.frames ready for pairwise comparisons
  df <- left_join(df_S, df_B[,colnames(df_B)!="scenario"], by = c("simulation_index"="simulation_index",
                                                                 "age_group"="age_group",
                                                                 "run_nr"="run_nr")) %>%
    rename(
      value_S = value.x,
      value_B = value.y
    )
  
  # Age stratified outcomes
  df_age <- df %>%
    group_by(scenario, simulation_index, age_group) %>%
    summarise(
      outcome_tot = median(value_B-value_S),
      outcome_per = median( (value_B-value_S)/value_B ),
      .groups = "drop"
    )
  
  # summed across all age-groups
  df_all <- df %>%
    filter(age_group == "All") %>%
    group_by(scenario, simulation_index) %>%
    summarise(
      outcome_tot = median(value_B-value_S),
      outcome_per = median( (value_B-value_S)/value_B ),
      .groups = "drop"
    )
  
  df_age_summary <- df_age %>%
    group_by(scenario, age_group) %>%
    summarise(
      median_tot = median(outcome_tot),
      median_per = median(outcome_per),
      upr95_per = quantile(outcome_per, 0.975),
      lwr95_per = quantile(outcome_per, 0.025),
      upr50_per = quantile(outcome_per, 0.75),
      lwr50_per = quantile(outcome_per, 0.25),
      .groups = "drop"
    )
  
  df_all_summary <- df_all %>%
    group_by(scenario) %>%
    summarise(
      median_tot = median(outcome_tot),
      median_per = median(outcome_per),
      upr95_per = quantile(outcome_per, 0.975),
      lwr95_per = quantile(outcome_per, 0.025),
      upr50_per = quantile(outcome_per, 0.75),
      lwr50_per = quantile(outcome_per, 0.25),
      .groups = "drop"
    )
  
  # Add labels for all the scenarios using look up tables
  df_all <- left_join(df_all, lookup_table, by = c("scenario" = "scenario"))
  df_age <- left_join(df_age, lookup_table, by = c("scenario" = "scenario"))
  df_all_summary <- left_join(df_all_summary, lookup_table, by = c("scenario" = "scenario"))
  df_age_summary <- left_join(df_age_summary, lookup_table, by = c("scenario" = "scenario"))  
  
  if (target == "infection") {
    df_inf_all <- df_all
    df_inf_age <- df_age
    df_inf_all_summary <- df_all_summary
    df_inf_age_summary <- df_age_summary
  } else if (target == "disease") {
    df_dis_all <- df_all
    df_dis_age <- df_age
    df_dis_all_summary <- df_all_summary
    df_dis_age_summary <- df_age_summary
  } else if (target == "admission") {
    df_adm_all <- df_all
    df_adm_age <- df_age
    df_adm_all_summary <- df_all_summary
    df_adm_age_summary <- df_age_summary
  } else if (target == "admission_best_flucan") {
    df_adm_best_flucan_all <- df_all
    df_adm_best_flucan_age <- df_age
    df_adm_best_flucan_all_summary <- df_all_summary
    df_adm_best_flucan_age_summary <- df_age_summary
  } else if (target == "admission_worst_flucan") {
    df_adm_worst_flucan_all <- df_all
    df_adm_worst_flucan_age <- df_age
    df_adm_worst_flucan_all_summary <- df_all_summary
    df_adm_worst_flucan_age_summary <- df_age_summary
  }
  
  cat("Finished loading", target, "data.\n\n")
}

rm(df, df_age, df_age_summary, df_all, df_all_summary, df_B, df_S, target)

##############################################################################################################

##############################################################################################################

# Figure 2
# Percentage infection reductions as a function of VE, per coverage percentage

fig_aggregated_results <- ggplot(df_inf_all[df_inf_all$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),], 
                                 aes(x=scen_effect, y= outcome_per, color=scen_coverage))+
  geom_point(position = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.75), alpha=0.15, size = 0.8, shape=16)+
  geom_pointrange(
    data = df_inf_all_summary[df_inf_all_summary$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),],
    aes(x=scen_effect, y=median_per, ymin=lwr50_per, ymax=upr50_per, group=scen_coverage),
    position = position_dodge(width = 0.75), inherit.aes = FALSE, colour = "black", linewidth = 0.65 , fatten = 2
    )+
  facet_wrap(~scen_age, ncol =1, labeller = as_labeller(c(
    "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
    "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
  )))+
  scale_color_brewer("LAIV\ncoverage",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(alpha = 1, size = 2)
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
      )+
  labs(
    x = "LAIV effectiveness",
    y = "Proportion of infections prevented"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    ),
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    panel.spacing = unit(0.35, "cm"),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

ggsave(filename = file.path(output_path,"perc_inf_red_for_ve_cov.png"), plot = fig_aggregated_results, height = 7.2, width = 6.3, units = "in", device = "png", dpi = 300)

##############################################################################################################

# Figure severe years

################################################################################
## Hospitalisations prevented by source
## 60% LAIV coverage only
################################################################################

# Combine the three sources into one dataframe
df_adm_sources <- bind_rows(
  df_adm_all %>%
    mutate(source = "AIHW"),
  
  df_adm_best_flucan_all %>%
    mutate(source = "FluCAN - low"),
  
  df_adm_worst_flucan_all %>%
    mutate(source = "FluCAN - high")
) %>%
  mutate(
    source = factor(
      source,
      levels = c("AIHW", "FluCAN - high", "FluCAN - low")
    )
  ) %>%
  filter(
    scen_age %in% c(
      "LAIV in 5-12 year olds",
      "LAIV in 5-18 year olds"
    ),
    scen_coverage == "60%"
  )

# Calculate summary statistics for number of hospitalisations prevented
df_adm_sources_summary <- df_adm_sources %>%
  group_by(scen_age, scen_effect, source) %>%
  summarise(
    median_tot = median(outcome_tot),
    lwr50_tot = quantile(outcome_tot, 0.25),
    upr50_tot = quantile(outcome_tot, 0.75),
    .groups = "drop"
  )

# Produce figure
fig_aggregated_admissions <- ggplot(
  df_adm_sources,
  aes(
    x = scen_effect,
    y = outcome_tot,
    color = source
  )
) +
  
  # Individual model particles
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.35,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.15,
    size = 0.8,
    shape = 16
  ) +
  
  # Median and 50% interval
  geom_pointrange(
    data = df_adm_sources_summary,
    aes(
      x = scen_effect,
      y = median_tot,
      ymin = lwr50_tot,
      ymax = upr50_tot,
      group = source
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.65,
    fatten = 2
  ) +
  
  # Separate panels for the two LAIV age groups
  facet_wrap(
    ~scen_age,
    ncol = 1,
    labeller = as_labeller(c(
      "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
      "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
    ))
  ) +
  
  # Colours now indicate source rather than coverage
  scale_color_brewer(
    "Hospitalisation data source",
    palette = "Dark2",
    labels = c(
      "AIHW" = "AIHW",
      "FluCAN - high" = "FluCAN adjusted – high",
      "FluCAN - low" = "FluCAN adjusted – low"
    ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  # Number of hospitalisations prevented
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02)),
    labels = scales::label_comma()
  ) +
  
  labs(
    x = "LAIV effectiveness",
    y = "Number of hospitalisations prevented per 100,000"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.3, "cm"),
    legend.margin = margin(t=2, r = 0, b = 0, l = 0),
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    legend.box.margin = margin(
      t = 5,
      r = 8,
      b = 5,
      l = 8
    ),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    panel.spacing = unit(0.35, "cm"),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# Save figure
ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_source_60coverage.png"
  ),
  plot = fig_aggregated_admissions,
  width = 6.3,
  height = 7.2,
  units = "in",
  dpi = 600
)

################################################################################
## Age-stratified hospitalisations prevented by source
## Central VE, 60% LAIV coverage
################################################################################

# Combine the three sources into one dataframe
df_adm_age_sources <- bind_rows(
  df_adm_age %>%
    mutate(source = "AIHW"),
  
  df_adm_best_flucan_age %>%
    mutate(source = "FluCAN - low"),
  
  df_adm_worst_flucan_age %>%
    mutate(source = "FluCAN - high")
) %>%
  mutate(
    source = factor(
      source,
      levels = c("AIHW", "FluCAN - high", "FluCAN - low")
    )
  ) %>%
  filter(
    scen_age %in% c(
      "LAIV in 5-12 year olds",
      "LAIV in 5-18 year olds"
    ),
    scen_coverage == "60%",
    scen_effect == "Central"
  )


# Make sure age groups appear in the correct order
df_adm_age_sources <- df_adm_age_sources %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c(
        "<1",
        "1-4",
        "5-11",
        "12-17",
        "18-64",
        "65-79",
        "80+",
        "All"
      )
    )
  )

# Calculate summary statistics
df_adm_age_sources_summary <- df_adm_age_sources %>%
  group_by(scen_age, age_group, source) %>%
  summarise(
    median_tot = median(outcome_tot),
    lwr50_tot = quantile(outcome_tot, 0.25),
    upr50_tot = quantile(outcome_tot, 0.75),
    .groups = "drop"
  )

# Produce figure
fig_age_stratified_admissions <- ggplot(
  df_adm_age_sources,
  aes(
    x = age_group,
    y = outcome_tot,
    color = source
  )
) +
  
  # Individual model particles
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.35,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.15,
    size = 0.8,
    shape = 16
  ) +
  
  # Median and 50% interval
  geom_pointrange(
    data = df_adm_age_sources_summary,
    aes(
      x = age_group,
      y = median_tot,
      ymin = lwr50_tot,
      ymax = upr50_tot,
      group = source
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.65,
    fatten = 2
  ) +
  geom_vline(
    xintercept = which(levels(factor(df_adm_age_sources$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  # Separate panels for the two LAIV age ranges
  facet_wrap(
    ~scen_age,
    ncol = 1,
    labeller = as_labeller(c(
      "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
      "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
    ))
  ) +
  
  # Colours indicate hospitalisation data source
  scale_color_brewer(
    "Hospitalisation data source",
    palette = "Dark2",
    labels = c(
      "AIHW" = "AIHW",
      "FluCAN - high" = "FluCAN adjusted – high",
      "FluCAN - low" = "FluCAN adjusted – low"
    ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  # Number of hospitalisations prevented
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02)),
    labels = scales::label_comma()
  ) +
  
  labs(
    x = "Age group (years)",
    y = "Number of hospitalisations prevented per 100,000"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_text(size = 11),
    axis.text.x = element_text(
      size = 10,
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.3, "cm"),
    legend.margin = margin(t=2, r = 0, b = 0, l = 0),
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    legend.box.margin = margin(
      t = 5,
      r = 8,
      b = 5,
      l = 8
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    panel.spacing = unit(0.35, "cm"),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# Save figure
ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_age_source_60coverage.png"
  ),
  plot = fig_age_stratified_admissions,
  width = 6.3,
  height = 7.2,
  units = "in",
  dpi = 600
)

##############################################################################################################

################################################################################
## Age-stratified proportion hospitalisations prevented by source
## Central VE, 60% LAIV coverage
################################################################################

# Combine the three sources into one dataframe
df_adm_age_sources <- bind_rows(
  df_adm_age %>%
    mutate(source = "AIHW"),
  
  df_adm_best_flucan_age %>%
    mutate(source = "FluCAN - low"),
  
  df_adm_worst_flucan_age %>%
    mutate(source = "FluCAN - high")
) %>%
  mutate(
    source = factor(
      source,
      levels = c("AIHW", "FluCAN - high", "FluCAN - low")
    )
  ) %>%
  filter(
    scen_age %in% c(
      "LAIV in 5-12 year olds",
      "LAIV in 5-18 year olds"
    ),
    scen_coverage == "60%",
    scen_effect == "Central"
  )


# Make sure age groups appear in the correct order
df_adm_age_sources <- df_adm_age_sources %>%
  mutate(
    age_group = factor(
      age_group,
      levels = c(
        "<1",
        "1-4",
        "5-11",
        "12-17",
        "18-64",
        "65-79",
        "80+", 
        "All"
      )
    )
  )


# Calculate summary statistics
df_adm_age_sources_summary <- df_adm_age_sources %>%
  group_by(scen_age, age_group, source) %>%
  summarise(
    median_tot = median(outcome_per),
    lwr50_tot = quantile(outcome_per, 0.25),
    upr50_tot = quantile(outcome_per, 0.75),
    .groups = "drop"
  )


# Produce figure
fig_age_stratified_admissions_perc <- ggplot(
  df_adm_age_sources,
  aes(
    x = age_group,
    y = outcome_per,
    color = source
  )
) +
  
  # Individual model particles
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.35,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.15,
    size = 0.8,
    shape = 16
  ) +
  
  # Median and 50% interval
  geom_pointrange(
    data = df_adm_age_sources_summary,
    aes(
      x = age_group,
      y = median_tot,
      ymin = lwr50_tot,
      ymax = upr50_tot,
      group = source
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.65,
    fatten = 2
  ) +
  geom_vline(
    xintercept = which(levels(factor(df_adm_age_sources$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  
  # Separate panels for the two LAIV age ranges
  facet_wrap(
    ~scen_age,
    ncol = 1,
    labeller = as_labeller(c(
      "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
      "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
    ))
  ) +
  
  # Colours indicate hospitalisation data source
  scale_color_brewer(
    "Hospitalisation data source",
    palette = "Dark2",
    labels = c(
      "AIHW" = "AIHW",
      "FluCAN - high" = "FluCAN adjusted – high",
      "FluCAN - low" = "FluCAN adjusted – low"
    ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  # Number of hospitalisations prevented
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02)),
    labels = scales::label_comma()
  ) +
  
  labs(
    x = "Age group (years)",
    y = "Proportion of hospitalisations prevented per 100,000"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_text(size = 11),
    axis.text.x = element_text(
      size = 10,
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.3, "cm"),
    legend.margin = margin(t=2, r = 0, b = 0, l = 0),
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    legend.box.margin = margin(
      t = 5,
      r = 8,
      b = 5,
      l = 8
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    panel.spacing = unit(0.35, "cm"),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# Save figure
ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_age_source_60coverage_perc.png"
  ),
  plot = fig_age_stratified_admissions_perc,
  width = 6.3,
  height = 7.2,
  units = "in",
  dpi = 600
)

##############################################################################################################

# Figure 3 Panel A
# Figure Infection reduction as a function of age group, per ve 

fig_age_results <- ggplot(df_inf_age[df_inf_age$scen_effect=="Central" &df_inf_age$scen_age %in% c("LAIV in 5-18 year olds"),], aes(x=age_group, y= outcome_per, color=scen_coverage))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_summary[df_inf_age_summary$scen_effect=="Central" &df_inf_age_summary$scen_age %in% c("LAIV in 5-18 year olds"),],
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=scen_coverage),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("LAIV coverage\nin 5 to <18 years",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(
    legend.background= element_rect(color="black"),
    strip.background = element_rect(color="black", fill="white"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = text_size_other),
    legend.text = element_text(size = text_size_other),
    
    axis.title.x = element_text(size = text_size_axis_title),
    axis.title.y = element_text(size = text_size_axis_title),
    
    axis.text.x = element_text(
      size = text_size_axis_text,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = text_size_axis_text,
      colour = "black"
    )
  )

# fig_age_results
# 
# ggsave(file.path(output_path, paste("inf_red_per_age_for_ve.png", sep = "")), height = 10, width = 14, dpi = 300)


##############################################################################################################
# Figure 3 panel B
# Negative particle plot - Panel instances of negatives

lookup_worst_reductions_idx <- df_inf_all %>%
  ungroup() %>%
  filter(scenario %in% c("Central_cover60_LAIV_5-18yo")) %>%
  arrange(outcome_per) %>%
  select(simulation_index, outcome_per) %>%
  filter(outcome_per < 0) %>%
  slice_head(n = 4) %>%
  pull(simulation_index) %>%
  .[-3] # We only want 3 particles, but the third one is too similar to one of the first two, so I'm filtering out the third one, and keeping numbers 1, 2, and 4

df_worst_particles <- particles %>%
  filter(simulation_index %in% lookup_worst_reductions_idx)%>%
  select(-c(frac_S1, half_life, phi_S1, phi_V0, phi_V1)) %>%
  crossing(t = t_vals) %>%
  mutate(
    beta = sin(2 * pi * (t - shift) / 365),
    simulation_index = factor(simulation_index)
  )

df_for_neg_path <- file.path(ROOT, paste("R outputs with waning", vaccination_term), "dat_uom_infection_with_waning")

df_worst_reductions <- open_dataset(df_for_neg_path) %>%
  filter(
    scenario %in% c("Status_quo", "Central_cover60_LAIV_5-18yo"), 
    simulation_index %in% lookup_worst_reductions_idx,
    run_nr == 0
  ) %>%
  select(-c(target, setting, team, run_nr)) %>%
  group_by(scenario, simulation_index, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect() %>%
  mutate(
    simulation_index = factor(
      simulation_index,
      levels = lookup_worst_reductions_idx
    )
  )

p1 <- ggplot(
  df_worst_reductions, 
  aes(
    x = horizon, 
    y = value, 
    colour = simulation_index,
    linetype = scenario, 
    group = interaction(simulation_index, scenario)
  )
) +
  geom_line(
    alpha = 1,
    linewidth = 1.2
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(
    name = NULL, breaks = c("Status_quo", "Central_cover60_LAIV_5-18yo"),
    values = c("Status_quo" = "solid", "Central_cover60_LAIV_5-18yo" = "dashed"),
    labels = c("Status_quo" = "Baseline", "Central_cover60_LAIV_5-18yo" = "Central 5 to <18 years")
  ) +
  guides(
    colour = "none",
    linetype = guide_legend(
      order = 1,
      override.aes = list(colour = "black")
    )
  ) +
  labs(x = NULL, y = "Infection incidence") +
  coord_cartesian(xlim = c(80, 300)) +
  theme_bw(base_size = text_size_other)+
  theme(
    # legend.position = c(0.98, 0.98),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = "black"
    ),
    
    legend.key.width = grid::unit(1.4, "cm"),
    legend.title = element_text(size = text_size_other),
    legend.text = element_text(size = text_size_other),
    
    strip.background = element_rect(color = "black", fill = "white"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = text_size_axis_text, colour = "black"),
    axis.title.y = element_text(size = text_size_axis_title),
    
    plot.margin = margin(8, 8, 2, 8)
  )

p2 <- ggplot(
  df_worst_particles, 
  aes(
    x = t, 
    y = beta, 
    colour = simulation_index,
    group = simulation_index
    )
) +
  geom_line(
    alpha = 1, 
    linetype = "dotdash",
    linewidth = 1.2
  ) + 
  scale_color_brewer(palette="Dark2")+
  guides(colour = "none") +
  labs(x = "Day", y = expression(beta(t))) +
  coord_cartesian(xlim = c(80, 300)) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = text_size_axis_text, colour = "black"),
        axis.text.y = element_text(size = text_size_axis_text, colour = "black"),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        plot.margin = margin(2, 8, 8, 8)
  )

fig_worst <- 
  (p1 / p2) +
  plot_layout(heights = c(5, 2))

# fig_worst
# 
# ggsave(file.path(output_path, paste("negative_reduction_analysis_plot.png", sep = "")), height = 10, width = 14, dpi = 300)

################################################################################
# Figure 3 panel C
# Negative particle plot - Panel peak time vs shift peak plot
# Version with only Baseline and Central 5-18, 60% coverage

source_file_timing_path <- file.path(ROOT, paste("Forward projection", source_file_type, vaccination_term), "df_max_timing.parquet")

scenarios_to_plot <- c(
  "Central_cover60_LAIV_5-18yo", "Status quo"
)

df_inf_for_colour <- df_inf_all %>%
  filter(scenario %in% c("Central_cover60_LAIV_5-18yo")) %>%
  select(simulation_index, outcome_per)

df_max_timing <- open_dataset(source_file_timing_path) %>%
  filter(scenario %in% scenarios_to_plot) %>%
  collect() %>%
  left_join(beta_peak_timing, by = "simulation_index") %>%
  left_join(df_inf_for_colour, by = "simulation_index") %>%
  rename(beta_peak_day = t)

min_val <- min(df_max_timing$outcome_per, na.rm = TRUE)
max_val <- max(df_max_timing$outcome_per, na.rm = TRUE)

make_timing_plot <- function(data, scenario_name, plot_title, show_x = TRUE, show_y = TRUE) {
  
  ggplot(
    data %>%
      filter(scenario == scenario_name), aes(x = max_day, y = beta_peak_day, colour = outcome_per)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
    annotate(
      "text",
      x = (120 + 365) / 2,
      y = 365,
      label = plot_title,
      hjust = 0.5,
      vjust = -0.45,
      size = text_size_other / ggplot2::.pt,
      lineheight = 0.85
    ) +
    coord_fixed(xlim = c(120, 365),ylim = c(120, 365), expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = seq(120, 360, by = 60)) +
    scale_y_continuous(breaks = seq(120, 360, by = 60)) +
    scale_colour_gradientn(
      colours = c(
        "#67000D",
        "#CB181D",
        "#FC9272",
        "grey98",
        "#9ECAE1",
        "#3182BD",
        "#08519C"
      ),
      values = scales::rescale(
        c(
          min_val,
          min_val * 0.6,
          min_val * 0.2,
          0,
          0.10,
          0.30,
          max_val
        ),
        from = c(min_val, max_val)
      ),
      limits = c(min_val, max_val),
      breaks = c(0, 0.25, 0.50, 0.75),
      labels = function(x) {
        ifelse(abs(x) < 1e-10, "0", sprintf("%.2f", x))
      },
      name = "Proportion of infections prevented in\nCentral 5 to <18 years compared to Baseline",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(15, "cm"),
        barheight = grid::unit(0.55, "cm"),
        title.theme = element_text(size = text_size_other, margin = margin(b = 2)),
        label.theme = element_text(size = text_size_axis_text)
      )
    ) +
    labs(
      x = if (show_x) "Peak infection incidence day" else NULL,
      y = if (show_y) "Peak transmissibility day" else NULL
    ) +
    theme_bw(base_size = text_size_other) +
    theme(
      axis.title.x = if (show_x) {
        element_text(size = text_size_axis_title, margin = margin(t = 5))
      } else {
        element_blank()
      },
      
      axis.title.y = if (show_y) {
        element_text(size = text_size_axis_title, margin = margin(r = 5))
      } else {
        element_blank()
      },
      
      axis.text.x = if (show_x) {
        element_text(
          size = text_size_axis_text,
          colour = "black"
        )
      } else {
        element_blank()
      },
      
      axis.text.y = if (show_y) {
        element_text(
          size = text_size_axis_text,
          colour = "black"
        )
      } else {
        element_blank()
      },
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.35
      ),
      
      panel.border = element_rect(
        colour = "grey35",
        linewidth = 0.7
      ),
      
      plot.margin = margin(
        t = 18,
        r = 12,
        b = 0,
        l = 6
      )
    )
}

timing_plot_layout <- tribble(
  ~scenario_name,                          ~title,                       ~row, ~column,
  "Status quo",                            "Baseline",                   1,       1,
  "Central_cover60_LAIV_5-18yo",           "Central 5 to <18 years",        1,      2
)

plots <- lapply(
  seq_len(nrow(timing_plot_layout)),
  function(i) {
    
    make_timing_plot(
      data = df_max_timing,
      scenario_name = timing_plot_layout$scenario_name[i],
      plot_title = timing_plot_layout$title[i],
      
      # Show x-axis only on the bottom row
      show_x = timing_plot_layout$row[i] == 1,
      
      # Show y-axis on the left column and on the baseline plot
      show_y = timing_plot_layout$column[i] == 1 
    )
  }
)

names(plots) <- timing_plot_layout$scenario_name

fig_timing <-
  patchwork::wrap_plots(
    plots[["Status quo"]],
    patchwork::plot_spacer(),
    plots[["Central_cover60_LAIV_5-18yo"]],
    nrow = 1,
    widths = c(1, 0.035, 1)
  ) +
  patchwork::plot_layout(
    guides = "collect"
    # widths = c(1,1)
  ) +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.box.background = element_rect(
      colour = "grey35",
      fill = NA,
      linewidth = 0.7
    ),
    legend.box.margin = margin(
      t = 0,
      r = 3,
      b = 3,
      l = 3
    ),
    legend.margin = margin(
      t = 0,
      r = 2,
      b = 2,
      l = 2
    ),
    legend.title = element_text(size = text_size_other, hjust = 0.5, margin = margin(b = 2)),
    legend.text = element_text(size = text_size_axis_text)
  )

# fig_timing
# 
# ggsave(file.path(output_path, paste("Diagonal_peak_plot_single_scenario.png", sep = "")), height = 7, width = 11, dpi = 300)


################################################################################

################################################################################

# Figure 3 - Combined panels

# Function for making a left-aligned panel label
make_panel_label <- function(label, parse = FALSE) {
  
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = label,
      parse = parse,
      hjust = 0,
      vjust = 0.5,
      size = text_size_panel_tag / ggplot2::.pt
      # fontface = "bold"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 2,
        l = 0
      )
    )
}

label_A <- make_panel_label(
  "bold('A: Proportion of infections prevented per age group')",
  parse = TRUE
)

label_B <- make_panel_label(
  "bold('B: Infection incidence and the corresponding ') * bold(beta(t))",
  parse = TRUE
)

label_C <- make_panel_label(
  "bold('C: Peak infection day vs. peak transmissibility day')",
  parse = TRUE
)

# Complete Panel A:
# label directly above the age-group figure
panel_A <-
  patchwork::wrap_plots(
    label_A,
    patchwork::wrap_elements(full = fig_age_results),
    ncol = 1,
    heights = c(0.06, 1)
  )

# Complete Panel B:
# label directly above fig_worst
panel_B <-
  patchwork::wrap_plots(
    label_B,
    patchwork::wrap_elements(full = fig_worst),
    ncol = 1,
    heights = c(0.06, 1)
  )

# Complete Panel C:
# label directly above fig_timing
panel_C <-
  patchwork::wrap_plots(
    label_C,
    patchwork::wrap_elements(full = fig_timing),
    ncol = 1,
    heights = c(0.06, 1)
  )

bottom_row <-
  patchwork::wrap_plots(
    panel_B,
    patchwork::plot_spacer(),
    panel_C,
    nrow = 1,
    widths = c(1, 0.08, 1)
  )

final_figure <-
  patchwork::wrap_plots(
    panel_A,
    patchwork::plot_spacer(),
    bottom_row,
    ncol = 1,
    heights = c(1, 0.08, 1)
  ) +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 10
      )
    )
  )

# final_figure

ggsave(
  filename = file.path(
    output_path,
    "combined_figures_A_B_C.png"
  ),
  plot = final_figure,
  width = 18,
  height = 13,
  dpi = 300,
  bg = "white"
)

################################################################################

################################################################################
# Figure 5
# Hospitalisation comparisons for good and bad years



################################################################################
# Figure 6
# Baseline peak beta sensitivity figure

df_inf_age_peak_beta <- df_inf_age %>%
  left_join(lookup_table_peak_beta, by = c("simulation_index")) %>%
  mutate(peak_beta_category = factor(peak_beta_category, levels=c("< 180", "180-210", "210-240", ">= 240"))) %>%
  filter(
    scen_effect == "Central" &
      scen_coverage == "60%" &
      scen_age == "LAIV in 5-18 year olds"
  )

df_inf_age_peak_beta_summary <- df_inf_age_peak_beta %>%
  group_by(scenario, age_group, peak_beta_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_peak_beta, aes(x=age_group, y= outcome_per, color=peak_beta_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_peak_beta_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=peak_beta_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Day of\nmaximal beta",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("Peak_beta_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

# rm(df_inf_age_peak_beta, df_inf_age_peak_beta_summary, lookup_table_peak_beta)

################################################################################

################################################################################

################################################################################
# Figure 6
# Waning sensitivity figure

df_inf_age_waning <- df_inf_age %>%
  left_join(lookup_table_half_life, by = c("simulation_index")) %>%
  mutate(half_life_category = factor(half_life_category, levels=c("<= 2y", "2y - 3y", "3y - 4y", "> 4y"))) %>%
  filter(
    scen_effect == "Central" &
    scen_coverage == "60%" &
      scen_age == "LAIV in 5-18 year olds"
    )

df_inf_age_waning_summary <- df_inf_age_waning %>%
  group_by(scenario, age_group, half_life_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_waning, aes(x=age_group, y= outcome_per, color=half_life_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_waning_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=half_life_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Waning\nhalf-life",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
        )

ggsave(file.path(output_path, paste("waning_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(df_inf_age_waning, df_inf_age_waning_summary, lookup_table_half_life)

################################################################################

################################################################################
# Figure 6
# Baseline AR sensitivity figure

df_inf_age_AR <- df_inf_age %>%
  left_join(lookup_table_AR, by = c("simulation_index")) %>%
  mutate(AR_category = factor(AR_category, levels=c("10-25", "25-35", "35-45"))) %>%
  filter(
    scen_effect == "Central" &
      scen_coverage == "60%" &
      scen_age == "LAIV in 5-18 year olds"
  )

df_inf_age_AR_summary <- df_inf_age_AR %>%
  group_by(scenario, age_group, AR_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_AR, aes(x=age_group, y= outcome_per, color=AR_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_AR_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=AR_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Baseline\nattack rate\n% of population",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("AR_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(df_inf_age_AR, df_inf_age_AR_summary, lookup_table_AR)

################################################################################

################################################################################
# Figure 6
# Baseline peak timing sensitivity figure

source_file_timing_path <- file.path(ROOT, paste("Forward projection", source_file_type, vaccination_term), "df_max_timing.parquet")

scenarios_to_plot <- c(
  "Status quo"
)

df_max_timing <- open_dataset(source_file_timing_path) %>%
  filter(scenario %in% scenarios_to_plot) %>%
  collect() %>%
  group_by(scenario, simulation_index) %>%
  summarise(max_day = median(max_day), .groups = "drop")

quantiles <- quantile(df_max_timing$max_day, probs = c(0.33, 0.66))

p33 <- quantiles[1]
p66 <- quantiles[2]

df_max_timing <- df_max_timing %>%
  mutate(peak_time_category = case_when(
    max_day < p33 ~ "Early",
    max_day >= p33 & max_day < p66 ~ "Mid-season",
    max_day >= p66 ~ "Late"
  )) %>%
  select(simulation_index, peak_time_category)


df_inf_age_peak_time <- df_inf_age %>%
  left_join(df_max_timing, by = c("simulation_index")) %>%
  mutate(peak_time_category = factor(peak_time_category, levels=c("Early", "Mid-season", "Late"))) %>%
  filter(
    scen_effect == "Central" &
      scen_coverage == "60%" &
      scen_age == "LAIV in 5-18 year olds"
  )

df_inf_age_peak_time_summary <- df_inf_age_peak_time %>%
  group_by(scenario, age_group, peak_time_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_peak_time, aes(x=age_group, y= outcome_per, color=peak_time_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_peak_time_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=peak_time_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Baseline\npeak time",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("peak_time_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(df_inf_age_peak_time, df_inf_age_peak_time_summary, source_file_timing_path, scenarios_to_plot, quantiles, p33, p66, df_max_timing)

################################################################################

################################################################################
# Figure 6
# Seasonality (beta_1 only) sensitivity figure

df_inf_age_seasonality <- df_inf_age %>%
  left_join(lookup_table_seasonality, by = c("simulation_index")) %>%
  mutate(beta_1_category = factor(beta_1_category, levels=c("weak", "medium", "strong"))) %>%
  mutate(beta0_x_beta1_category = factor(beta0_x_beta1_category, levels=c("weak", "medium", "strong"))) %>%
  filter(
    scen_effect == "Central" &
      scen_coverage == "60%" &
      scen_age == "LAIV in 5-18 year olds"
  )

df_inf_age_seasonality_beta1_summary <- df_inf_age_seasonality %>%
  group_by(scenario, age_group, beta_1_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_seasonality, aes(x=age_group, y= outcome_per, color=beta_1_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_seasonality_beta1_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=beta_1_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Seasonality\nparameter beta_1",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("seasonality_beta1_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

df_inf_age_seasonality_beta0_x_beta1_summary <- df_inf_age_seasonality %>%
  group_by(scenario, age_group, beta0_x_beta1_category) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

ggplot(df_inf_age_seasonality, aes(x=age_group, y= outcome_per, color=beta0_x_beta1_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_seasonality_beta0_x_beta1_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=beta0_x_beta1_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("Seasonality parameter\n beta_0 x beta_1",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(
                         alpha = 1,
                         size = 2
                       )
                     ))+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("seasonality_beta0_x_beta1_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(df_inf_age_seasonality, df_inf_age_seasonality_beta1_summary, df_inf_age_seasonality_beta0_x_beta1_summary, lookup_table_seasonality)

################################################################################

################################################################################
# Figure 7
# Effect of timing of vaccination

source_file_path <- file.path(ROOT, paste("R outputs", source_file_type, "term2 vaccination"), "infection_ar_ds")

df_inf_term2 <- open_dataset(source_file_path) %>%
  collect()

df_term2_all_ages <- df_inf_term2 %>%
  group_by(scenario, simulation_index, run_nr) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(age_group = "All")

df_inf_term2 <- dplyr::bind_rows(df_inf_term2, df_term2_all_ages)
rm(df_term2_all_ages)

df_inf_term2$age_group <- factor(df_inf_term2$age_group, levels = age_levels)

## Summarising
cat("Loading infection data for term 2 vaccination run.\n")

# Seperate dataframes for baseline vs scenarios
df_B <- df_inf_term2[df_inf_term2$scenario=="Status_quo",]
df_S <- df_inf_term2[df_inf_term2$scenario!="Status_quo",]  

# Join data.frames ready for pairwise comparisons
df_inf_term2 <- left_join(df_S, df_B[colnames(df_B)!="scenario"], by = c("simulation_index"="simulation_index",
                                                                         "age_group"="age_group",
                                                                         "run_nr"="run_nr")) %>%
  rename(
    value_S = value.x,
    value_B = value.y
  )

# Age stratified outcomes
df_inf_age_term2 <- df_inf_term2 %>%
  group_by(scenario, simulation_index, age_group) %>%
  summarise(
    outcome_tot = median(value_B-value_S),
    outcome_per = median( (value_B-value_S)/value_B ),
    .groups = "drop"
  )

# summed across all age-groups
df_inf_all_term2 <- df_inf_term2 %>%
  filter(age_group == "All") %>%
  group_by(scenario, simulation_index) %>%
  summarise(
    outcome_tot = median(value_B-value_S),
    outcome_per = median( (value_B-value_S)/value_B ),
    .groups = "drop"
  )

df_inf_age_summary_term2 <- df_inf_age_term2 %>%
  group_by(scenario, age_group) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

df_inf_all_summary_term2 <- df_inf_all_term2 %>%
  group_by(scenario) %>%
  summarise(
    median_tot = median(outcome_tot),
    median_per = median(outcome_per),
    upr95_per = quantile(outcome_per, 0.975),
    lwr95_per = quantile(outcome_per, 0.025),
    upr50_per = quantile(outcome_per, 0.75),
    lwr50_per = quantile(outcome_per, 0.25),
    .groups = "drop"
  )

# Add labels for all the scenarios using look up tables
df_inf_all_term2 <- left_join(df_inf_all_term2, lookup_table, by = c("scenario" = "scenario"))
df_inf_age_term2 <- left_join(df_inf_age_term2, lookup_table, by = c("scenario" = "scenario"))
df_inf_all_summary_term2 <- left_join(df_inf_all_summary_term2, lookup_table, by = c("scenario" = "scenario"))
df_inf_age_summary_term2 <- left_join(df_inf_age_summary_term2, lookup_table, by = c("scenario" = "scenario"))  

cat("Finished loading infection data for term 2 vaccination run.\n\n")

rm(df_B, df_S, lookup_table)

df_inf_term1 <- df_inf
df_inf_age_term1 <- df_inf_age
df_inf_age_summary_term1 <- df_inf_age_summary
df_inf_all_term1 <- df_inf_all
df_inf_all_summary_term1 <- df_inf_all_summary

rm(df_inf, df_inf_age, df_inf_age_summary, df_inf_all, df_inf_all_summary)

df_inf_term1_adj <- df_inf_age_term1
df_inf_term2_adj <- df_inf_age_term2
df_inf_term1_adj$label <- "Term 1"
df_inf_term2_adj$label <- "Term 2"
df_timing <- rbind(df_inf_term1_adj, df_inf_term2_adj)

sum_term1 <- df_inf_age_summary_term1
sum_term2 <- df_inf_age_summary_term2
sum_term1$label <- "Term 1"
sum_term2$label <- "Term 2"
df_timing_summary <- rbind(sum_term1, sum_term2)

ggplot(df_timing[df_timing$scen_effect=="Central" &df_timing$scen_age %in% c("LAIV in 5-18 year olds") &df_timing$scen_coverage=="60%",], aes(x=age_group, y= outcome_per, color=label))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_timing_summary[df_timing_summary$scen_effect=="Central" &df_timing_summary$scen_age %in% c("LAIV in 5-18 year olds")&df_timing_summary$scen_coverage=="60%",,],
                  aes(x=age_group, y=median_per,ymin=upr50_per,ymax=lwr50_per,  group=label),
                  size=1,
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,
                  color="black")+
  scale_color_brewer("LAIV\ntiming",palette="Dark2")+
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0)),
    breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_vline(
    xintercept = which(levels(factor(df_timing$age_group)) == "All") - 0.5,
    linewidth = 1,
    colour = "black"
  ) +
  theme_bw(base_size = text_size_other)+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = text_size_other),
        legend.text = element_text(size = text_size_other),
        
        axis.title.x = element_text(size = text_size_axis_title),
        axis.title.y = element_text(size = text_size_axis_title),
        
        axis.text.x = element_text(
          size = text_size_axis_text,
          colour = "black"
        ),
        axis.text.y = element_text(
          size = text_size_axis_text,
          colour = "black"
        )
  )

ggsave(file.path(output_path, paste("perc_inf_averted_per_timing_LAIV.png", sep = "")), height = 10, width = 14, dpi = 300)

################################################################################
# Supplements figures
################################################################################
# Figure 4 - Supplement
# Negative particle plot - Panel peak time vs shift peak plot
# Version with all 60% scenarios

scenarios_to_plot <- c(
  "Status quo", 
  "Pessimistic_cover60_LAIV_5-12yo", "Pessimistic_cover60_LAIV_5-18yo",
  "Central_cover60_LAIV_5-12yo", "Central_cover60_LAIV_5-18yo",
  "Optimistic_cover60_LAIV_5-12yo", "Optimistic_cover60_LAIV_5-18yo"
)

df_max_timing <- open_dataset(source_file_timing_path) %>%
  filter(scenario %in% scenarios_to_plot) %>%
  collect() %>%
  left_join(beta_peak_timing, by = "simulation_index") %>%
  left_join(df_inf_for_colour, by = "simulation_index") %>%
  rename(beta_peak_day = t)

min_val <- min(df_max_timing$outcome_per, na.rm = TRUE)
max_val <- max(df_max_timing$outcome_per, na.rm = TRUE)

make_timing_plot <- function(data, scenario_name, plot_title, show_x = TRUE, show_y = TRUE) {
  
  ggplot(
    data %>%
      filter(scenario == scenario_name), aes(x = max_day, y = beta_peak_day, colour = outcome_per)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
    coord_fixed(xlim = c(120, 365),ylim = c(120, 365), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 360, by = 60)) +
    scale_y_continuous(breaks = seq(0, 360, by = 60)) +
    scale_colour_gradientn(
      colours = c(
        "#67000D",
        "#CB181D",
        "#FC9272",
        "grey98",
        "#9ECAE1",
        "#3182BD",
        "#08519C"
      ),
      values = scales::rescale(
        c(
          min_val,
          min_val * 0.6,
          min_val * 0.2,
          0,
          0.10,
          0.30,
          max_val
        ),
        from = c(min_val, max_val)
      ),
      limits = c(min_val, max_val),
      name = "Reduction \nproportion"
    ) +
    labs(
      title = plot_title,
      x = if (show_x) "Peak infection day" else NULL,
      y = if (show_y) "Peak transmission day" else NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, margin = margin(b = 8)),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11, colour = "black"),
      axis.text.x = if (show_x) {
        element_text(size = 11, colour = "black")
      } else {
        element_blank()
      },
      axis.title.x = if (show_x) {
        element_text(size = 13)
      } else {
        element_blank()
      },
      axis.text.y = if (show_y) {
        element_text(size = 11, colour = "black")
      } else {
        element_blank()
      },
      axis.title.y = if (show_y) {
        element_text(size = 13)
      } else {
        element_blank()
      },
      
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.35
      ),
      
      panel.border = element_rect(
        colour = "grey35",
        linewidth = 0.7
      ),
      
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_layout <- tribble(
  ~scenario_name,                          ~title,                       ~row, ~column,
  "Status quo",                            "Baseline",                   1,       2,
  "Pessimistic_cover60_LAIV_5-12yo",       "Pessimistic\n5–12 years",    2,       1,
  "Central_cover60_LAIV_5-12yo",           "Central\n5–12 years",        2,       2,
  "Optimistic_cover60_LAIV_5-12yo",        "Optimistic\n5–12 years",     2,       3,
  "Pessimistic_cover60_LAIV_5-18yo",       "Pessimistic\n5–18 years",    3,       1,
  "Central_cover60_LAIV_5-18yo",           "Central\n5–18 years",        3,       2,
  "Optimistic_cover60_LAIV_5-18yo",        "Optimistic\n5–18 years",     3,       3
)

plots <- lapply(
  seq_len(nrow(plot_layout)),
  function(i) {
    
    make_timing_plot(
      data = df_max_timing,
      scenario_name = plot_layout$scenario_name[i],
      plot_title = plot_layout$title[i],
      
      # Show x-axis only on the bottom row
      show_x = plot_layout$row[i] == 3,
      
      # Show y-axis on the left column and on the baseline plot
      show_y = plot_layout$column[i] == 1 |
        plot_layout$scenario_name[i] == "Status quo"
    )
  }
)

names(plots) <- plot_layout$scenario_name

fig_timing <-
  (
    plot_spacer() |
      plots[["Status quo"]] |
      plot_spacer()
  ) /
  (
    plots[["Pessimistic_cover60_LAIV_5-12yo"]] |
      plots[["Central_cover60_LAIV_5-12yo"]] |
      plots[["Optimistic_cover60_LAIV_5-12yo"]]
  ) /
  (
    plots[["Pessimistic_cover60_LAIV_5-18yo"]] |
      plots[["Central_cover60_LAIV_5-18yo"]] |
      plots[["Optimistic_cover60_LAIV_5-18yo"]]
  ) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.margin = margin(15, 15, 15, 15)
    )
  ) &
  theme(
    legend.position = "bottom",
    
    legend.box.background = element_rect(
      colour = "grey35",
      fill = NA,
      linewidth = 0.7
    ),
    
    legend.box.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 10
    ),
    
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.width = grid::unit(4, "cm")
  )

fig_timing

ggsave(file.path(output_path, paste("Diagonal_peak_plot_all_cover60.png", sep = "")), height = 11, width = 11, dpi = 300)


################################################################################
# Older versions of figures that may still be useful
################################################################################
# Figure 4 with A and B vertically stacked

combined_figure_vertical <- patchwork::wrap_plots(
  patchwork::wrap_elements(full = fig_worst),
  patchwork::wrap_elements(full = fig_timing),
  ncol = 1,
  heights = c(1, 0.9)
) +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  )

# combined_figure_vertical

ggsave(file.path(output_path, paste("combined_vertical.png", sep = "")),  plot = combined_figure_vertical, height = 17, width = 12, dpi = 300)


############

################################################################################
## Percentage reduction in hospitalisations by source
## 60% LAIV coverage only
################################################################################

# Combine the three sources into one dataframe
df_adm_sources_per <- bind_rows(
  df_adm_all %>%
    mutate(source = "AIHW"),
  
  df_adm_best_flucan_all %>%
    mutate(source = "FluCAN - best"),
  
  df_adm_worst_flucan_all %>%
    mutate(source = "FluCAN - worst")
) %>%
  filter(
    scen_age %in% c(
      "LAIV in 5-12 year olds",
      "LAIV in 5-18 year olds"
    ),
    scen_coverage == "60%"
  )


# Calculate summary statistics
df_adm_sources_per_summary <- df_adm_sources_per %>%
  group_by(scen_age, scen_effect, source) %>%
  summarise(
    median_per = median(outcome_per),
    lwr50_per = quantile(outcome_per, 0.25),
    upr50_per = quantile(outcome_per, 0.75),
    .groups = "drop"
  )


# Produce figure
fig_aggregated_admissions_per <- ggplot(
  df_adm_sources_per,
  aes(
    x = scen_effect,
    y = outcome_per,
    color = source
  )
) +
  
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.35,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.15,
    size = 0.8,
    shape = 16
  ) +
  
  geom_pointrange(
    data = df_adm_sources_per_summary,
    aes(
      x = scen_effect,
      y = median_per,
      ymin = lwr50_per,
      ymax = upr50_per,
      group = source
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.65,
    fatten = 2
  ) +
  
  facet_wrap(
    ~scen_age,
    ncol = 1,
    labeller = as_labeller(c(
      "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
      "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
    ))
  ) +
  
  scale_color_brewer(
    "Hospitalisation\ndata source",
    palette = "Dark2",
    guide = guide_legend(
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult = c(0.02, 0.01)),
    breaks = seq(-0.2, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  
  labs(
    x = "LAIV effectiveness",
    y = "Reduction in hospitalisations"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    ),
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    
    panel.spacing = unit(0.35, "cm"),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )


ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_reduction_percentage_by_source_60coverage.png"
  ),
  plot = fig_aggregated_admissions_per,
  height = 7.2,
  width = 6.3,
  units = "in",
  device = "png",
  dpi = 300
)

