rm(list = ls())

source_file_type <- "with waning"
vaccination_term <- "term1 vaccination"

target_list <- c("infection", "disease", "admission")

library(ggplot2)
library(patchwork)
library(arrow)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)

source_file_formatted <- gsub("\\s+", "_", source_file_type)

age_levels <- c("<1", "1-4", "5-11", "12-17", "18-64", "65-79", "80+")

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
  
  df$age_group <- factor(df$age_group, levels = age_levels)
  
  if (target == "infection") {
    df_inf <- df
  } else if (target == "disease") {
    df_dis <- df
  } else if (target == "admission") {
    df_adm <- df
  }
}

particles_file_path <- file.path(ROOT, paste("Forward projection with waning", vaccination_term), "particles.csv")
particles <- read.csv(particles_file_path) %>%
  select(-X) %>%
  rownames_to_column(var = "simulation_index") %>%
  mutate(simulation_index = as.integer(simulation_index)-1)

rm(df, particles_file_path, target, source_file_path)

#################################################################################################################
# This section to be commented out if term 2 vaccination is not also being loaded in

source_file_path <- file.path(ROOT, paste("R outputs", source_file_type, "term2 vaccination"), "infection_ar_ds")

df_inf_term2 <- open_dataset(source_file_path) %>%
  collect()

df_inf_term2$age_group <- factor(df_inf_term2$age_group, levels = age_levels)

rm(source_file_path)

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
  group_by(scenario, simulation_index, run_nr) %>%
  summarise(AR = sum(value), .groups = "drop") %>%
  group_by(scenario, simulation_index) %>%
  summarise(AR = median(AR), .groups = "drop")  %>%
  filter(scenario == "Status_quo")%>%
  mutate(AR_category = case_when(
    AR >= 100000 & AR < 250000 ~ "10-25",
    AR >= 250000 & AR < 350000 ~ "25-35",
    AR >= 350000               ~ "35-45" 
  )) %>%
  select(simulation_index, AR_category)

lookup_table_seasonality <- particles %>%
  select(simulation_index, beta_0, beta_1) %>%
  mutate(
    beta_0_x_beta_1 = beta_0 * beta_1,
    beta_1_category = case_when(
      beta_1 < 0.133 ~ "weak",
      beta_1 >= 0.133 & beta_1 < 0.216 ~ "medium",
      beta_1 >= 0.216 ~ "strong"
    ),
    beta0_x_beta1_category = case_when(
      beta_0_x_beta_1 < 0.00667 ~ "weak",
      beta_0_x_beta_1 >= 0.00667 & beta_0_x_beta_1 < 0.0113 ~ "medium",
      beta_0_x_beta_1 >= 0.0113 ~ "strong"
    )
    ) %>%
  select(-c(beta_0, beta_1, beta_0_x_beta_1))

t_vals <- seq(0, 364)
beta_peak_timing <- particles %>%
  crossing(t = t_vals) %>%
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

rm(scen_age, scen_coverage, scen_effect, scenarios_considered)

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
  }
  
  # Seperate dataframes for baseline vs scenarios
  df_B <- df[df$scenario=="Status_quo",]
  df_S <- df[df$scenario!="Status_quo",]  
  
  # Join data.frames ready for pairwise comparisons
  df <- left_join(df_S, df_B[colnames(df_B)!="scenario"], by = c("simulation_index"="simulation_index",
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
    group_by(scenario, simulation_index, run_nr) %>%
    summarise(
      value_B = sum(value_B),
      value_S = sum(value_S),
      .groups = "drop"
    ) %>%
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
  }
  
  cat("Finished loading", target, "data.\n\n")
}

rm(df, df_age, df_age_summary, df_all, df_all_summary, df_B, df_S, target)

##############################################################################################################

################################################################################
# Comment out this section if term 2 vaccination is not being run
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
  group_by(scenario, simulation_index, run_nr) %>%
  summarise(
    value_B = sum(value_B),
    value_S = sum(value_S),
    .groups = "drop"
  ) %>%
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

##############################################################################################################

# Figure 2
# Percentage infection reductions as a function of VE, per coverage percentage

ggplot(df_inf_all[df_inf_all$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),], aes(x=scen_effect, y= outcome_per, color=scen_coverage))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(
    data = df_inf_all_summary[df_inf_all_summary$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),],
    aes(x=scen_effect, y=median_per, ymin=lwr50_per, ymax=upr50_per, group=scen_coverage),
    position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
    shape=95, size=1, color="black")+
  facet_wrap(.~scen_age, nrow=3)+
  scale_color_brewer("LAIV\ncoverage",palette="Dark2",
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
  xlab("LAIV effectiveness")+
  ylab("Proportion infections prevented")+
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.y = element_blank())

ggsave(file.path(output_path, paste("perc_inf_red_for_ve_cov.png", sep = "")), height = 10, width = 14, dpi = 300)

##############################################################################################################
# Figure 3
# Figure Infection reduction as a function of age group, per ve 

ggplot(df_inf_age[df_inf_age$scen_effect=="Central" &df_inf_age$scen_age %in% c("LAIV in 5-18 year olds"),], aes(x=age_group, y= outcome_per, color=scen_coverage))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, shape=16)+
  geom_pointrange(data= df_inf_age_summary[df_inf_age_summary$scen_effect=="Central" &df_inf_age_summary$scen_age %in% c("LAIV in 5-18 year olds"),],
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=scen_coverage),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  shape=95,size=1, color="black")+
  scale_color_brewer("LAIV\ncoverage",palette="Dark2",
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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave(file.path(output_path, paste("inf_red_per_age_for_ve.png", sep = "")), height = 10, width = 14, dpi = 300)


##############################################################################################################
# Figure 4 A
# Negative particle plot - Panel instances of negatives

lookup_worst_reductions_idx <- df_inf_all %>%
  ungroup() %>%
  filter(scenario %in% c("Central_cover60_LAIV_5-18yo")) %>%
  arrange(outcome_per) %>%
  mutate(improvement_rank = row_number()) %>%
  select(simulation_index, outcome_per, improvement_rank) %>%
  filter(
    outcome_per < 0
  ) %>%
  slice_head(n = 4) %>%
  pull(simulation_index) %>%
  .[-3] # We only want 3 particles, but the third one is too similar to one of the first two, so I'm filtering out the third one, and keeping numbers 1, 2, and 4

df_worst_particles <- particles[particles$simulation_index %in% lookup_worst_reductions_idx,] %>%
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
  collect()

p1 <- ggplot(
  df_worst_reductions, 
  aes(
    x = horizon, 
    y = value, 
    colour = factor(simulation_index),
    linetype = scenario, 
    group = interaction(simulation_index, scenario)
  )
) +
  geom_line(
    # colour = "black",
    alpha = 1,
    linewidth = 1.2
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(
    name = NULL, breaks = c("Status_quo", "Central_cover60_LAIV_5-18yo"),
    values = c("Status_quo" = "solid", "Central_cover60_LAIV_5-18yo" = "dashed"),
    labels = c("Status_quo" = "Baseline", "Central_cover60_LAIV_5-18yo" = "With LAIV")
  ) +
  guides(colour = "none")+
  labs(x = NULL, y = "Infections") +
  coord_cartesian(xlim = c(80, 300)) +
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )

p2 <- ggplot(
  df_worst_particles, 
  aes(
    x = t, 
    y = beta, 
    colour = factor(simulation_index),
    group = simulation_index)
) +
  geom_line(
    # colour = "black",
    alpha = 1, 
    linetype = "dotdash",
    linewidth = 1.2
  ) + 
  scale_color_brewer(palette="Dark2")+
  guides(colour = "none") +
  labs(x = "Day", y = expression(beta(t))) +
  coord_cartesian(xlim = c(80, 300)) +
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()
  )

(p1 / p2) +
  plot_layout(
    heights = c(5, 2),
    guides = "collect"
  ) &
  theme(legend.position = "bottom",
        legend.key.width=grid::unit(1.5, "cm"))


ggsave(file.path(output_path, paste("negative_reduction_analysis_plot.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(lookup_worst_reductions_idx, df_worst_particles, df_for_neg_path, df_worst_reductions, p1, p2)

################################################################################
# Figure 4B
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
    coord_fixed(xlim = c(0, 365),ylim = c(0, 365), expand = FALSE) +
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
      name = "Infection reduction \nproportion"
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
  "Central_cover60_LAIV_5-18yo",           "Central\n5–18 years",        1,       1
)

plots <- lapply(
  seq_len(nrow(plot_layout)),
  function(i) {
    
    make_timing_plot(
      data = df_max_timing,
      scenario_name = plot_layout$scenario_name[i],
      plot_title = plot_layout$title[i],
      
      # Show x-axis only on the bottom row
      show_x = plot_layout$row[i] == 1,
      
      # Show y-axis on the left column and on the baseline plot
      show_y = plot_layout$column[i] == 1 
    )
  }
)

names(plots) <- plot_layout$scenario_name

fig_timing <-
  (
    plots[["Central_cover60_LAIV_5-18yo"]] |
      plots[["Status quo"]] 
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

ggsave(file.path(output_path, paste("Diagonal_peak_plot_single_scenario.png", sep = "")), height = 11, width = 11, dpi = 300)

################################################################################
# Figure 4B
# Negative particle plot - Panel peak time vs shift peak plot
# Version with all 60% scenarios

source_file_timing_path <- file.path(ROOT, paste("Forward projection", source_file_type, vaccination_term), "df_max_timing.parquet")

scenarios_to_plot <- c(
  "Status quo", 
  "Pessimistic_cover60_LAIV_5-12yo", "Pessimistic_cover60_LAIV_5-18yo",
  "Central_cover60_LAIV_5-12yo", "Central_cover60_LAIV_5-18yo",
  "Optimistic_cover60_LAIV_5-12yo", "Optimistic_cover60_LAIV_5-18yo"
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
    coord_fixed(xlim = c(0, 365),ylim = c(0, 365), expand = FALSE) +
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
      name = "Infection reduction \nproportion"
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

rm(df_inf_for_colour, min_val, max_val, make_timing_plot, fig_timing, plots, plot_layout, beta_peak_timing)

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave(file.path(output_path, paste("seasonality_beta0_x_beta1_sensitivity_of_infections_central_scenario.png", sep = "")), height = 10, width = 14, dpi = 300)

rm(df_inf_age_seasonality, df_inf_age_seasonality_beta1_summary, df_inf_age_seasonality_beta0_x_beta1_summary, lookup_table_seasonality)

################################################################################

################################################################################
# Figure 7
# Effect of timing of vaccination

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
  theme_bw()+
  theme(legend.background= element_rect(color="black"),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

ggsave(file.path(output_path, paste("perc_inf_averted_per_timing_LAIV.png", sep = "")), height = 10, width = 14, dpi = 300)






