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
age_levels_no_all <- c("<1", "1-4", "5-11", "12-17", "18-64", "65-79", "80+")

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

common_theme <- theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    
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
    
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    
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
                          "Optimistic_cover80_LAIV_5-18yo"
                          # "Central_cover40_LAIV_2-5yo"
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
              "LAIV in 5-18 year olds"
              # "LAIV in 2-5 year olds"
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
                 "Optimistic"
                 # "Central"
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
                   "80%"
                   # "40%"
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

lookup_table$scen_age <- factor(
  lookup_table$scen_age,
  levels = c(
    "Baseline",
    "LAIV in 5-12 year olds",
    "LAIV in 5-18 year olds"
    # "LAIV in 2-5 year olds"
  )
)

lookup_table$scen_coverage <- factor(
  lookup_table$scen_coverage,
  levels = c("Baseline", "20%", "40%", "60%", "80%")
)

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
  
  df_age$age_group <- factor(df_age$age_group, levels=age_levels)
  df_age_summary$age_group <- factor(df_age_summary$age_group, levels=age_levels)
  
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

##############################################################################################################
# FIGURE INFECTION REDUCTION AS A FUNCTION OF VE, PER COVERAGE PERCENTAGE
##############################################################################################################

fig_aggregated_results <- ggplot(df_inf_all[df_inf_all$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),], 
                                 aes(x=scen_effect, y= outcome_per, color=scen_coverage))+
  geom_point(position = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.75), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(
    data = df_inf_all_summary[df_inf_all_summary$scen_age %in% c("LAIV in 5-12 year olds","LAIV in 5-18 year olds"),],
    aes(x=scen_effect, y=median_per, ymin=lwr50_per, ymax=upr50_per, group=scen_coverage),
    position = position_dodge(width = 0.75), inherit.aes = FALSE, colour = "black", linewidth = 0.65 , fatten = 2
    )+
  facet_wrap(~scen_age, ncol =1, labeller = as_labeller(c(
    "LAIV in 5-12 year olds" = "LAIV 5 to <12 years",
    "LAIV in 5-18 year olds" = "LAIV 5 to <18 years"
  )))+
  scale_color_brewer("LAIV coverage",palette="Dark2",
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
  
  # theme_bw(base_size = 10) +
  
  common_theme +
  theme(
    # legend.position = "right",
    # legend.direction = "vertical",
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    )
  ) 

ggsave(filename = file.path(output_path,"perc_inf_red_for_ve_cov.png"), plot = fig_aggregated_results, height = 8.38, width = 7.33, units = "in", device = "png", dpi = 300)


##############################################################################################################
##############################################################################################################

##############################################################################################################
# FIGURE EFFECT OF ROLLOUT TIMING
##############################################################################################################

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

fig_vacc_term_effect <- ggplot(df_timing[df_timing$scen_effect=="Central" &df_timing$scen_age %in% c("LAIV in 5-18 year olds") &df_timing$scen_coverage=="60%",], aes(x=age_group, y= outcome_per, color=label))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(data= df_timing_summary[df_timing_summary$scen_effect=="Central" &df_timing_summary$scen_age %in% c("LAIV in 5-18 year olds")&df_timing_summary$scen_coverage=="60%",,],
                  aes(x=age_group, y=median_per,ymin=upr50_per,ymax=lwr50_per,  group=label),
                  # size=1,
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  # shape=95,
                  color="black")+
  
  scale_color_brewer(
    "Rollout timing",
    palette = "Dark2",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
    )+
  
  labs(
    x = "Age group (years)",
    y = "Proportion of infections prevented"
  )+
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  )+
  geom_vline(
    xintercept = which(levels(factor(df_timing$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
    common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )

ggsave(file.path(output_path, paste("perc_inf_averted_per_timing_LAIV.png", sep = "")), fig_vacc_term_effect, height = 5.0, width = 7.33, dpi = 300)

####################################################
# END PERCENTAGE INFECTIONS AVERTED PER TERM LAIV
####################################################


###############################################################################
# TIMING OF BASELINE EPIDEMIC TRAJECTORIES RELATIVE TO THE VACCINATION WINDOW
###############################################################################

source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning") 
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_daily_inf_filtered_sims <- ds_daily_inf %>%
  filter(
    scenario %in% c("Status_quo"),
    # simulation_index < 250,
    run_nr == 0
  ) %>%
  select(c(scenario, simulation_index, age_group, horizon, run_nr, value)) %>%
  group_by(simulation_index, run_nr, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect() 

fig_daily_inf_filtered_sims <- ggplot(
  df_daily_inf_filtered_sims,
  aes(
    x = horizon,
    y = value,
    group = interaction(simulation_index, run_nr)
  )
) +
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3,
    colour = "royalblue"
  ) +
  
  scale_y_continuous(
    labels = function(x) x / 1e3
  ) +
  
  annotate(
    "rect",
    xmin = 110,
    xmax = 170,
    ymin = -Inf,
    ymax = Inf,
    fill = "#1B9E77",
    alpha = 0.3
  ) +
  
  annotate(
    "rect",
    xmin = 59,
    xmax = 100,
    ymin = -Inf,
    ymax = Inf,
    fill = "#D95F02",
    alpha = 0.3
  ) +
  
  labs(
    x = "Day",
    y = "Daily infection incidence ('000s)",
  ) +
  
  common_theme +
  
  theme(
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    output_path,
    "daily_inf_inc_vs_vacc_window.png"
  ),
  plot = fig_daily_inf_filtered_sims,
  height = 4.5,
  width = 7.33,
  units = "in",
  device = "png",
  dpi = 300
)

###################################################################################
# END TIMING OF BASELINE EPIDEMIC TRAJECTORIES RELATIVE TO THE VACCINATION WINDOW
###################################################################################

##############################################################################################################
# FIGURES SEVERE YEARS
##############################################################################################################

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
    alpha = 0.2,
    size = 1.0,
    shape = 16
  ) +
  
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
    "Hospitalisation burden",
    palette = "Dark2",
    labels = c(
      "AIHW" = "Primary analysis (AIHW 2023)",
      "FluCAN - high" = "High burden",
      "FluCAN - low" = "Low burden"
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
    expand = expansion(mult=c(0.02, 0.01)),
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
  
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    )
  ) 
  
# Save figure
ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_source_60coverage.png"
  ),
  plot = fig_aggregated_admissions,
  width = 7.33,
  height = 8.38,
  units = "in",
  device = "png",
  dpi = 300
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
    alpha = 0.2,
    size = 1.0,
    shape = 16
  ) +
  
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
    linewidth = 0.8,
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
  
  scale_color_brewer(
    "Hospitalisation burden",
    palette = "Dark2",
    labels = c(
      "AIHW" = "Primary analysis (AIHW 2023)",
      "FluCAN - high" = "High burden",
      "FluCAN - low" = "Low burden"
    ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  # Number of hospitalisations prevented
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.01)),
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
  
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

# Save figure
ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_age_source_60coverage.png"
  ),
  plot = fig_age_stratified_admissions,
  width = 7.33,
  height = 8.38,
  units = "in",
  device = "png",
  dpi = 300
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

df_adm_age_sources_summary <- df_adm_age_sources %>%
  group_by(scen_age, age_group, source) %>%
  summarise(
    median_tot = median(outcome_per),
    lwr50_tot = quantile(outcome_per, 0.25),
    upr50_tot = quantile(outcome_per, 0.75),
    .groups = "drop"
  )

fig_age_stratified_admissions_perc <- ggplot(
  df_adm_age_sources,
  aes(
    x = age_group,
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
    alpha = 0.2,
    size = 1.0,
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
    linewidth = 0.8,
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
  
  scale_color_brewer(
    "Hospitalisation burden",
    palette = "Dark2",
    labels = c(
      "AIHW" = "Primary analysis",
      "FluCAN - high" = "High burden",
      "FluCAN - low" = "Low burden"
    ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  # Number of hospitalisations prevented
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.01)),
    labels = scales::label_comma()
  ) +
  
  labs(
    x = "Age group (years)",
    y = "Proportion of hospitalisations prevented"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
    common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    output_path,
    "hospitalisations_prevented_by_age_source_60coverage_perc.png"
  ),
  plot = fig_age_stratified_admissions_perc,
  width = 7.33,
  height = 8.38,
  units = "in",
  device = "png",
  dpi = 300
)

##############################################################################################################
##############################################################################################################

##############################################################################################################
# FIGURE AGE-STRATIEFIED REDUCTIONS AND NEGATIVE RESULTS
##############################################################################################################

##############################################################################################################
# PANEL A
##############################################################################################################

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
    legend.title = element_text(size = text_size_legend_title),
    legend.text = element_text(size = text_size_legend_text),
    
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

ggsave(file.path(output_path, paste("inf_red_per_age_for_ve.png", sep = "")), fig_age_results, height = 10, width = 14, dpi = 300)


##############################################################################################################
# PANEL B
##############################################################################################################

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
      fill = alpha("white", 0.5),
      colour = "black"
    ),
    
    legend.key.width = grid::unit(1.4, "cm"),
    legend.title = element_text(size = text_size_legend_title),
    legend.text = element_text(size = text_size_legend_text+1),
    
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
# PANEL C
################################################################################
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
      name = "Proportion of infections prevented in\nCentral 5 to <18 years",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(15, "cm"),
        barheight = grid::unit(0.55, "cm"),
        title.theme = element_text(size = text_size_legend_title, margin = margin(b = 2)),
        label.theme = element_text(size = text_size_legend_text)
      )
    ) +
    labs(
      x = if (show_x) "Peak infection day" else NULL,
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
    legend.title = element_text(size = text_size_legend_title, hjust = 0.5, margin = margin(b = 2)),
    legend.text = element_text(size = text_size_legend_text)
  )

# fig_timing
# 
# ggsave(file.path(output_path, paste("Diagonal_peak_plot_single_scenario.png", sep = "")), height = 7, width = 11, dpi = 300)

################################################################################
# COMBINED PANELS
################################################################################

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
# SUPPLEMENTS
################################################################################

##############################################################################################################
# BASELINE INFECTION ATTACK RATE
##############################################################################################################
df_baseline_only <- df_inf %>%
  filter(scenario == "Status_quo")

df_baseline_only$age_group <- factor(df_baseline_only$age_group, levels = age_levels)

fig_baseline_inf_ar <- ggplot(df_baseline_only, aes(x = age_group, y = value)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), colour = "royalblue", alpha=0.2, size = 1.0, shape=16)+
  
  scale_x_discrete(
    labels = c(
      "<1"     = "<1",
      "1-4"   = "1–<5",
      "5-11"  = "5–<12",
      "12-17" = "12–<18",
      "18-64" = "18–<65",
      "65-79" = "65–<80",
      "80+"   = "80+",
      "All"   = "All"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult=c(0.02, 0.01)),
    labels = function(x) x / 1e5
  )+
  
  labs(
    x = "Age group (years)",
    y = "Infection attack rate ('00000s)"
  )+
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  )+
  geom_vline(
    xintercept = which(levels(factor(df_baseline_only$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  common_theme +
  theme(
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )

ggsave(file.path(output_path, paste("baseline_infection_ar.png", sep = "")), fig_baseline_inf_ar, height = 5.0, width = 7.33, dpi = 300)

##############################################################################################################

##############################################################################################################
# BASELINE INFECTION ATTACK RATE PER 100,000 PER AGE GROUP
##############################################################################################################
df_baseline_only <- df_inf %>%
  filter(scenario == "Status_quo") %>%
  left_join(ages, by = "age_group") %>%
  mutate(ipht = 100000 * value / pop_size)

df_baseline_only$age_group <- factor(df_baseline_only$age_group, levels = age_levels)

fig_baseline_inf_ar_per_ht <- ggplot(df_baseline_only, aes(x = age_group, y = ipht)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), colour = "royalblue", alpha=0.2, size = 1.0, shape=16)+
  
  scale_x_discrete(
    labels = c(
      "<1"     = "<1",
      "1-4"   = "1–<5",
      "5-11"  = "5–<12",
      "12-17" = "12–<18",
      "18-64" = "18–<65",
      "65-79" = "65–<80",
      "80+"   = "80+",
      "All"   = "All"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult=c(0.02, 0.01)),
    labels = function(x) x / 1e4
  )+
  
  labs(
    x = "Age group (years)",
    y = "Infection attack rate per 100,000 ('0000s)"
  )+
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  )+
  geom_vline(
    xintercept = which(levels(factor(df_baseline_only$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  common_theme +
  theme(
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )

ggsave(file.path(output_path, paste("baseline_infection_ar_per_100k.png", sep = "")), fig_baseline_inf_ar_per_ht, height = 5.0, width = 7.33, dpi = 300)

##############################################################################################################

##############################################################################################################
# SENTITIVITY ANALYSIS FIGURES
##############################################################################################################

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

fig_max_beta_sensitivity <- ggplot(df_inf_age_peak_beta, aes(x=age_group, y= outcome_per, color=peak_beta_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  
  geom_pointrange(data= df_inf_age_peak_beta_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=peak_beta_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  # shape=95,
                  color="black")+
  
  scale_color_brewer(
    expression("Day of maximal " * beta),
    # "Day of maximal beta",
    palette="Dark2",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
      )
    )+
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    # breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
    labels = scales::label_number(accuracy = 0.1)) +
  
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
  
  geom_hline(yintercept = 0, linetype="dashed", linewidth = 0.45, colour="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +

  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("Peak_beta_sensitivity_of_infections_central_scenario.png", sep = "")), fig_max_beta_sensitivity, height = 5.0, width = 7.33, dpi = 300)

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

fig_waning_sensitivity <- ggplot(df_inf_age_waning, aes(x=age_group, y= outcome_per, color=half_life_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  
  geom_pointrange(data= df_inf_age_waning_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=half_life_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  colour = "black")+
  
  scale_color_brewer("Waning half-life",palette="Dark2",
                     guide = guide_legend(
                       direction = "horizontal",
                       title.position = "top",
                       title.hjust = 0.5,
                       override.aes = list(alpha = 1, size = 2)
                       )
                     )+
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1))+
  
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("waning_sensitivity_of_infections_central_scenario.png", sep = "")), fig_waning_sensitivity, height = 5.0, width = 7.33, dpi = 300)

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

fig_ar_sensitivity <- ggplot(df_inf_age_AR, aes(x=age_group, y= outcome_per, color=AR_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(data= df_inf_age_AR_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=AR_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  color="black")+
  
  scale_color_brewer("Baseline attack rate (% of population)",palette="Dark2",
                     guide = guide_legend(
                       direction = "horizontal",
                       title.position = "top",
                       title.hjust = 0.5,
                       override.aes = list(alpha = 1, size = 2)
                       )
                     )+
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)) +
  
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  
  geom_hline(yintercept = 0, linetype="dashed", linewidth = 0.45, color="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("AR_sensitivity_of_infections_central_scenario.png", sep = "")), fig_ar_sensitivity, height = 5.0, width = 7.33, dpi = 300)

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

fig_peak_timing_sensitivity <- ggplot(df_inf_age_peak_time, aes(x=age_group, y= outcome_per, color=peak_time_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(data= df_inf_age_peak_time_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=peak_time_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  color="black")+
  
  scale_color_brewer("Baseline peak time",palette="Dark2",
                     guide = guide_legend(
                       direction = "horizontal",
                       title.position = "top",
                       title.hjust = 0.5,
                       override.aes = list(alpha = 1, size = 2)
                       )
                     )+
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)) +
  
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("peak_time_sensitivity_of_infections_central_scenario.png", sep = "")), fig_peak_timing_sensitivity, height = 5.0, width = 7.33, dpi = 300)

# rm(df_inf_age_peak_time, df_inf_age_peak_time_summary, source_file_timing_path, scenarios_to_plot, quantiles, p33, p66, df_max_timing)

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

fig_beta_seasonality_sensitivity <- ggplot(df_inf_age_seasonality, aes(x=age_group, y= outcome_per, color=beta_1_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(data= df_inf_age_seasonality_beta1_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=beta_1_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  color="black")+
  
  scale_color_brewer(
    expression("Seasonality parameter " * beta[1]),
    palette="Dark2",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  )+
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)) +
  
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("seasonality_beta1_sensitivity_of_infections_central_scenario.png", sep = "")),fig_beta_seasonality_sensitivity, height = 5.0, width = 7.33, dpi = 300)

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

fig_beta_x_beta_seasonality_sensitivity <- ggplot(df_inf_age_seasonality, aes(x=age_group, y= outcome_per, color=beta0_x_beta1_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(data= df_inf_age_seasonality_beta0_x_beta1_summary,
                  aes(x=age_group, y=median_per, ymin=lwr50_per, ymax=upr50_per,  group=beta0_x_beta1_category),
                  position = position_jitterdodge(jitter.width = 0.0, jitter.height = 0),
                  inherit.aes = FALSE,
                  linewidth = 0.65,
                  fatten = 2,
                  color="black")+
  
  scale_color_brewer(
    expression("Seasonality parameter " * beta[0] %*% beta[1]),
    ,palette="Dark2",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  )+
  
  scale_x_discrete(
    labels = c(
      "1-4" = "1-<5",
      "5-11" = "5-<12",
      "12-17" = "12-<18",
      "18-64" = "18-<65",
      "65-79" = "65-<80"
    )
  ) +
  
  scale_y_continuous(
    limits = c(NA, 1),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks= seq(-0.2, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)) +
  
  xlab("Age-group")+
  ylab("Proportion of infections prevented")+
  
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  
  geom_vline(
    xintercept = which(levels(factor(df_inf_age$age_group)) == "All") - 0.5,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    
    axis.text.x = element_text(
      # size = 10,
      angle = 30,
      hjust = 1
    )
  )  

ggsave(file.path(output_path, paste("seasonality_beta0_x_beta1_sensitivity_of_infections_central_scenario.png", sep = "")), fig_beta_x_beta_seasonality_sensitivity, height = 5.0, width = 7.33, dpi = 300)

################################################################################
# STOCHASTICITY FIGURES
################################################################################

#########################################################################
# SELECTED PARTICLES - DAILY INFECTION INCIDENCE - AGE STRATIFIED
#########################################################################

scenario_levels <- c("Status_quo", "Central_cover60_LAIV_5-18yo", "Optimistic_cover60_LAIV_5-18yo")

source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning") 
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_daily_inf_selected_parts <- ds_daily_inf %>%
  filter(
    scenario %in% c("Status_quo", "Central_cover60_LAIV_5-18yo", "Optimistic_cover60_LAIV_5-18yo"),
    simulation_index %in% c(0, 7)
  ) %>%
  select(c(scenario, simulation_index, age_group, horizon, run_nr, value)) %>%
  collect() %>%
  mutate(
    age_group = factor(age_group, levels = age_levels_no_all),
    scenario = factor(scenario, levels = scenario_levels)
  )

fig_daily_inf_age_strat <- ggplot(
  df_daily_inf_selected_parts,
  aes(
    x = horizon,
    y = value,
    colour = factor(simulation_index),
    group = interaction(simulation_index, run_nr)
  )
) +
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3
  ) +
  
  facet_grid(
    rows = vars(age_group),
    cols = vars(scenario),
    scales = "free_y",
    labeller = labeller(
      scenario = c(
        "Status_quo" = "Baseline",
        "Central_cover60_LAIV_5-18yo" = "Central 5–18",
        "Optimistic_cover60_LAIV_5-18yo" = "Optimistic 5–18"
      ),
      age_group = c(
        "<1" = "<1",
        "1-4" = "1–<5",
        "5-11" = "5–<12",
        "12-17" = "12–<18",
        "18-64" = "18–<65",
        "65-79" = "65–<80",
        "80+" = "80+"
      )
    )
  ) +
  
  scale_x_continuous(
    breaks = c(0, 150, 300)
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "#1B9E77",
      "7" = "#D95F02"
    )
  ) +
  
  labs(
    x = "Day",
    y = "Daily infection incidence",
    # colour = NULL
  ) +
  
  common_theme +
  
  theme(
    # No x-axis labels or ticks
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    
    # Slightly tighter spacing because there are seven panels
    panel.spacing = unit(0.15, "cm"),
    
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    output_path,
    "selected_parts_stochasticity_age_stratified.png"
  ),
  plot = fig_daily_inf_age_strat,
  height = 8.9,
  width = 7.33,
  units = "in",
  device = "png",
  dpi = 300
)


#########################################################################
# SELECTED PARTICLES - DAILY INFECTION INCIDENCE - AGE AGGREGATED
#########################################################################

source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning") 
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_daily_inf_selected_parts <- ds_daily_inf %>%
  filter(
    simulation_index %in% c(0, 7, 56, 367)
  ) %>%
  select(scenario, simulation_index, age_group, horizon, run_nr, value) %>%
  group_by(scenario, simulation_index, horizon, run_nr) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect() %>%
  left_join(lookup_table, by = "scenario")

df_baseline <- df_daily_inf_selected_parts %>%
  filter(scenario == "Status_quo")

df_intervention <- df_daily_inf_selected_parts %>%
  filter(scenario != "Status_quo") %>%
  mutate(
    row_label = paste(
      scen_effect,
      ifelse(
        scen_age == "LAIV in 5-12 year olds",
        "5–12",
        "5–18"
      ),
      sep = "\n"
    ),
    row_label = factor(
      row_label,
      levels = c(
        "Pessimistic\n5–12",
        "Pessimistic\n5–18",
        "Central\n5–12",
        "Central\n5–18",
        "Optimistic\n5–12",
        "Optimistic\n5–18"
      )
    )
  )

fig_intervention <- ggplot(
  df_intervention,
  aes(
    x = horizon,
    y = value,
    colour = factor(simulation_index),
    group = interaction(simulation_index, run_nr)
  )
) +
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3
  ) +
  
  facet_grid(
    rows = vars(row_label),
    cols = vars(scen_coverage)
    ) +

  scale_x_continuous(
    breaks = c(0, 150, 300)
  ) +
  
  scale_y_continuous(
    labels = function(x) x / 1e3,
    breaks = c(0, 5000, 10000)
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "#1B9E77",
      "7" = "#D95F02",
      "56" = "#7570B3",
      "367" = "#E7298A"
    )
  ) +
  
  labs(
    x = "Day",
    y = "Daily infection incidence ('000s)"
  ) +
  
  common_theme +
  
  theme(
    panel.spacing = unit(0.15, "cm"),
    legend.position = "none"
    )
  # )

fig_baseline <- ggplot(
  df_baseline,
  aes(
    x = horizon,
    y = value,
    colour = factor(simulation_index),
    group = interaction(simulation_index, run_nr)
  )
) +
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3
  ) +
  
  scale_x_continuous(
    breaks = c(0, 150, 300)
  ) +
  
  scale_y_continuous(
    breaks = c(0, 5000, 10000),
    labels = function(x) x / 1e3
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "#1B9E77",
      "7" = "#D95F02",
      "56" = "#7570B3",
      "367" = "#E7298A"
    )
  ) +
  
  labs(
    title = "Baseline",
    x = NULL,
    y = NULL
  ) +
  
  common_theme +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

baseline_row <-
  plot_spacer() +
  fig_baseline +
  plot_spacer() +
  plot_layout(
    ncol = 3,
    widths = c(1.5, 1, 1.5)
  )

fig_daily_inf_all_scenarios <-
  baseline_row /
  fig_intervention +
  plot_layout(
    heights = c(1, 6)
  )

ggsave(
  filename = file.path(
    output_path,
    "daily_infections_all_scenarios.png"
  ),
  plot = fig_daily_inf_all_scenarios,
  width = 7.33,
  height = 8.9,
  units = "in",
  device = "png",
  dpi = 300
)


#########################################################################
# SELECTED PARTICLES - YEARLY ATTACK RATE PER SCENARIO - AGE STRATIFIED
#########################################################################

scenario_levels <- c(
  "Status_quo", 
  
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
  "Optimistic_cover80_LAIV_5-18yo"
)

df_ar_selected_parts_age_strat <- df_inf %>%
  filter(
    simulation_index %in% c(0, 7, 56, 367),
    age_group != "All"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  )

fig_inf_selected_parts_age_strat <- ggplot(
  df_ar_selected_parts_age_strat,
  aes(
    x = scenario,
    y = value,
    colour = factor(simulation_index),
    group = interaction(simulation_index, run_nr)
  )
) +
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3
  ) +
  geom_point(
    alpha = 0.2,
    size = 1.0,
    shape = 16
  ) +
  
  facet_grid(
    rows = vars(age_group),
    scales = "free_y",
    labeller = labeller(
      age_group = c(
        "<1" = "<1",
        "1-4" = "1–<5",
        "5-11" = "5–<12",
        "12-17" = "12–<18",
        "18-64" = "18–<65",
        "65-79" = "65–<80",
        "80+" = "80+"
      )
    )
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "#1B9E77",
      "7" = "#D95F02",
      "56" = "#7570B3",
      "367" = "#E7298A"
    )
  ) +
  
  labs(
    x = "Scenario",
    y = "Infection attack rate"
  ) +
  
  common_theme +
  
  theme(
    # No x-axis labels or ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Slightly tighter spacing because there are seven panels
    panel.spacing = unit(0.15, "cm"),
    
    legend.position = "none"
  )
  
  # scale_x_discrete(
  #   labels = c(
  #     "Status_quo" = "Baseline",
  #     
  #     "Pessimistic_cover20_LAIV_5-12yo" = "Pessimistic\n5–12, 20%",
  #     "Pessimistic_cover40_LAIV_5-12yo" = "Pessimistic\n5–12, 40%",
  #     "Pessimistic_cover60_LAIV_5-12yo" = "Pessimistic\n5–12, 60%",
  #     "Pessimistic_cover80_LAIV_5-12yo" = "Pessimistic\n5–12, 80%",
  #     "Pessimistic_cover20_LAIV_5-18yo" = "Pessimistic\n5–18, 20%",
  #     "Pessimistic_cover40_LAIV_5-18yo" = "Pessimistic\n5–18, 40%",
  #     "Pessimistic_cover60_LAIV_5-18yo" = "Pessimistic\n5–18, 60%",
  #     "Pessimistic_cover80_LAIV_5-18yo" = "Pessimistic\n5–18, 80%",
  #     
  #     "Central_cover20_LAIV_5-12yo" = "Central\n5–12, 20%",
  #     "Central_cover40_LAIV_5-12yo" = "Central\n5–12, 40%",
  #     "Central_cover60_LAIV_5-12yo" = "Central\n5–12, 60%",
  #     "Central_cover80_LAIV_5-12yo" = "Central\n5–12, 80%",
  #     "Central_cover20_LAIV_5-18yo" = "Central\n5–18, 20%",
  #     "Central_cover40_LAIV_5-18yo" = "Central\n5–18, 40%",
  #     "Central_cover60_LAIV_5-18yo" = "Central\n5–18, 60%",
  #     "Central_cover80_LAIV_5-18yo" = "Central\n5–18, 80%",
  #     
  #     "Optimistic_cover20_LAIV_5-12yo" = "Optimistic\n5–12, 20%",
  #     "Optimistic_cover40_LAIV_5-12yo" = "Optimistic\n5–12, 40%",
  #     "Optimistic_cover60_LAIV_5-12yo" = "Optimistic\n5–12, 60%",
  #     "Optimistic_cover80_LAIV_5-12yo" = "Optimistic\n5–12, 80%",
  #     "Optimistic_cover20_LAIV_5-18yo" = "Optimistic\n5–18, 20%",
  #     "Optimistic_cover40_LAIV_5-18yo" = "Optimistic\n5–18, 40%",
  #     "Optimistic_cover60_LAIV_5-18yo" = "Optimistic\n5–18, 60%",
  #     "Optimistic_cover80_LAIV_5-18yo" = "Optimistic\n5–18, 80%"
  #   )
  # ) +

ggsave(
  filename = file.path(
    output_path,
    "selected_parts_stochasticity_ar_per_scenario_age_strat.png"
  ),
  plot = fig_inf_selected_parts_age_strat,
  height = 8.9,
  width = 7.33,
  units = "in",
  device = "png",
  dpi = 300
)

#########################################################################
# SELECTED PARTICLES - YEARLY ATTACK RATE PER SCENARIO - AGE AGGREGATED
#########################################################################

scenario_levels <- c(
  "Status_quo", 
  
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
  "Optimistic_cover80_LAIV_5-18yo"
)

df_ar_selected_parts <- df_inf %>%
  filter(
    simulation_index %in% c(0, 7, 56, 367),
    age_group == "All"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  )

fig_inf_selected_parts <- ggplot(
  df_ar_selected_parts,
  aes(
    x = scenario,
    y = value,
    colour = factor(simulation_index),
    group = interaction(simulation_index, run_nr)
  )
) +
  
  scale_y_continuous(
    expand = expansion(mult=c(0.02, 0.01)),
    labels = function(x) x / 1e5
  )+
  
  geom_line(
    alpha = 0.2,
    linewidth = 0.3
  ) +
  geom_point(
    alpha = 0.2,
    size = 1.0,
    shape = 16
  ) +
  scale_colour_manual(
    values = c(
      "0" = "#1B9E77",
      "7" = "#D95F02",
      "56" = "#7570B3",
      "367" = "#E7298A"
    )
  ) +
  labs(
    x = "Scenario",
    y = "Infection attack rate ('00000s)"
  ) +
  # scale_x_discrete(
  #   labels = c(
  #     "Status_quo" = "Baseline",
  #     
  #     "Pessimistic_cover20_LAIV_5-12yo" = "Pessimistic\n5–12, 20%",
  #     "Pessimistic_cover40_LAIV_5-12yo" = "Pessimistic\n5–12, 40%",
  #     "Pessimistic_cover60_LAIV_5-12yo" = "Pessimistic\n5–12, 60%",
  #     "Pessimistic_cover80_LAIV_5-12yo" = "Pessimistic\n5–12, 80%",
  #     "Pessimistic_cover20_LAIV_5-18yo" = "Pessimistic\n5–18, 20%",
  #     "Pessimistic_cover40_LAIV_5-18yo" = "Pessimistic\n5–18, 40%",
  #     "Pessimistic_cover60_LAIV_5-18yo" = "Pessimistic\n5–18, 60%",
  #     "Pessimistic_cover80_LAIV_5-18yo" = "Pessimistic\n5–18, 80%",
  #     
  #     "Central_cover20_LAIV_5-12yo" = "Central\n5–12, 20%",
  #     "Central_cover40_LAIV_5-12yo" = "Central\n5–12, 40%",
  #     "Central_cover60_LAIV_5-12yo" = "Central\n5–12, 60%",
  #     "Central_cover80_LAIV_5-12yo" = "Central\n5–12, 80%",
  #     "Central_cover20_LAIV_5-18yo" = "Central\n5–18, 20%",
  #     "Central_cover40_LAIV_5-18yo" = "Central\n5–18, 40%",
  #     "Central_cover60_LAIV_5-18yo" = "Central\n5–18, 60%",
  #     "Central_cover80_LAIV_5-18yo" = "Central\n5–18, 80%",
  #     
  #     "Optimistic_cover20_LAIV_5-12yo" = "Optimistic\n5–12, 20%",
  #     "Optimistic_cover40_LAIV_5-12yo" = "Optimistic\n5–12, 40%",
  #     "Optimistic_cover60_LAIV_5-12yo" = "Optimistic\n5–12, 60%",
  #     "Optimistic_cover80_LAIV_5-12yo" = "Optimistic\n5–12, 80%",
  #     "Optimistic_cover20_LAIV_5-18yo" = "Optimistic\n5–18, 20%",
  #     "Optimistic_cover40_LAIV_5-18yo" = "Optimistic\n5–18, 40%",
  #     "Optimistic_cover60_LAIV_5-18yo" = "Optimistic\n5–18, 60%",
  #     "Optimistic_cover80_LAIV_5-18yo" = "Optimistic\n5–18, 80%"
  #   )
  # ) +
  
  common_theme +
  
  theme(
    # No x-axis labels or ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Slightly tighter spacing because there are seven panels
    panel.spacing = unit(0.15, "cm"),
    
    legend.position = "none"
  )
  
  # theme_bw() +
  # theme(
  #   axis.text.x = element_text(
  #     size = 7,
  #     angle = 30,
  #     hjust = 1,
  #     vjust = 1
  #   ),
  #   axis.text.y = element_text(size = 7),
  #   axis.title.y = element_text(size = 9),
  #   panel.grid.minor = element_blank(),
  #   legend.position = "none"
  # )

ggsave(
  filename = file.path(
    output_path,
    "selected_parts_stochastic_ar_per_scenario.png"
  ),
  plot = fig_inf_selected_parts,
  width = 7.33,
  height = 4.0,
  units = "in",
  device = "png",
  dpi = 300
)

################################################################################
# END STOCHASTICITY FIGURES
################################################################################

#########################################################################
# FIGURE - RANKING SCENARIOS BASED ON INFECTION REDUCTION
#########################################################################

scenario_levels <- c(
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
  "Optimistic_cover80_LAIV_5-18yo"
)

df_ar_all_parts <- df_inf %>%
  filter(
    age_group == "All"
  )

df_mean_ar <- df_ar_all_parts %>%
  group_by(
    scenario,
    simulation_index
  ) %>%
  summarise(
    mean_attack_rate = mean(value),
    .groups = "drop"
  )

df_baseline_ar <- df_mean_ar %>%
  filter(
    scenario == "Status_quo"
  ) %>%
  select(
    simulation_index,
    baseline_attack_rate = mean_attack_rate
  )

df_rank <- df_mean_ar %>%
  filter(
    scenario != "Status_quo"
  ) %>%
  left_join(
    df_baseline_ar,
    by = "simulation_index"
  ) %>%
  mutate(
    mean_ar_reduction =
      (baseline_attack_rate - mean_attack_rate) /
      baseline_attack_rate,
    
    scenario = factor(
      scenario,
      levels = scenario_levels
    )
  ) %>%
  group_by(
    scenario
  ) %>%
  mutate(
    rank = rank(
      -mean_ar_reduction,
      ties.method = "average"
    )
  ) %>%
  ungroup()

particle_order <- df_rank %>%
  filter(
    scenario == "Central_cover60_LAIV_5-18yo"
  ) %>%
  arrange(
    rank
  ) %>%
  pull(
    simulation_index
  )

df_rank <- df_rank %>%
  mutate(
    simulation_index = factor(
      simulation_index,
      levels = rev(particle_order)
    )
  )

fig_rank <- ggplot(
  df_rank,
  aes(
    x = scenario,
    y = simulation_index,
    fill = rank
  )
) +
  geom_tile() +
  
  scale_fill_distiller(
    palette = "RdYlBu",
    direction = -1
  ) +
  
  scale_x_discrete(
    labels = c(
      "Pessimistic_cover20_LAIV_5-12yo" = "Pessimistic\n5–12, 20%",
      "Pessimistic_cover40_LAIV_5-12yo" = "Pessimistic\n5–12, 40%",
      "Pessimistic_cover60_LAIV_5-12yo" = "Pessimistic\n5–12, 60%",
      "Pessimistic_cover80_LAIV_5-12yo" = "Pessimistic\n5–12, 80%",
      
      "Pessimistic_cover20_LAIV_5-18yo" = "Pessimistic\n5–18, 20%",
      "Pessimistic_cover40_LAIV_5-18yo" = "Pessimistic\n5–18, 40%",
      "Pessimistic_cover60_LAIV_5-18yo" = "Pessimistic\n5–18, 60%",
      "Pessimistic_cover80_LAIV_5-18yo" = "Pessimistic\n5–18, 80%",
      
      "Central_cover20_LAIV_5-12yo" = "Central\n5–12, 20%",
      "Central_cover40_LAIV_5-12yo" = "Central\n5–12, 40%",
      "Central_cover60_LAIV_5-12yo" = "Central\n5–12, 60%",
      "Central_cover80_LAIV_5-12yo" = "Central\n5–12, 80%",
      
      "Central_cover20_LAIV_5-18yo" = "Central\n5–18, 20%",
      "Central_cover40_LAIV_5-18yo" = "Central\n5–18, 40%",
      "Central_cover60_LAIV_5-18yo" = "Central\n5–18, 60%",
      "Central_cover80_LAIV_5-18yo" = "Central\n5–18, 80%",
      
      "Optimistic_cover20_LAIV_5-12yo" = "Optimistic\n5–12, 20%",
      "Optimistic_cover40_LAIV_5-12yo" = "Optimistic\n5–12, 40%",
      "Optimistic_cover60_LAIV_5-12yo" = "Optimistic\n5–12, 60%",
      "Optimistic_cover80_LAIV_5-12yo" = "Optimistic\n5–12, 80%",
      
      "Optimistic_cover20_LAIV_5-18yo" = "Optimistic\n5–18, 20%",
      "Optimistic_cover40_LAIV_5-18yo" = "Optimistic\n5–18, 40%",
      "Optimistic_cover60_LAIV_5-18yo" = "Optimistic\n5–18, 60%",
      "Optimistic_cover80_LAIV_5-18yo" = "Optimistic\n5–18, 80%"
    )
  ) +
  
  labs(
    x = "Scenario",
    y = "Particle",
    # fill = "Rank"
  ) +
  
  common_theme +
  
  theme(
    # No x and y -axes labels or ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # # Slightly tighter spacing because there are seven panels
    # panel.spacing = unit(0.15, "cm"),
    
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    output_path,
    "particle_rankings_ar_reduction_all_scenarios.png"
  ),
  plot = fig_rank,
  width = 7.33,
  height = 8.9,
  units = "in",
  device = "png",
  dpi = 300
)

################################################################################
# END RANKING SCENARIOS FIGURE
################################################################################

###############################################################################
# COMPARISON OF PEAK TIMING BETWEEN BASELINE, TERM 1 VACCINATION
# AND TERM 2 VACCINATION
###############################################################################

source_file_daily_inf_term1 <- file.path(ROOT, paste("R outputs", source_file_type, "term1 vaccination"), "dat_uom_infection_with_waning")
source_file_daily_inf_term2 <- file.path(ROOT, paste("R outputs", source_file_type, "term2 vaccination"), "dat_uom_infection_with_waning") 

ds_daily_inf_term1 <- open_dataset(source_file_daily_inf_term1)
ds_daily_inf_term2 <- open_dataset(source_file_daily_inf_term2)

df_daily_inf_filtered_sims_term1 <- ds_daily_inf_term1 %>%
  filter(
    scenario %in% c("Status_quo", "Central_cover60_LAIV_5-18yo"),
  ) %>%
  select(c(scenario, simulation_index, run_nr, age_group, horizon, value)) %>%
  group_by(scenario, simulation_index, run_nr, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect() %>%
  group_by(scenario, simulation_index, run_nr) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(peak_day = horizon) %>%
  mutate(
    type = case_when(
      scenario == "Central_cover60_LAIV_5-18yo"  ~ "term1",
      scenario == "Status_quo"  ~ "Baseline"
    )
  ) %>%
  select(type, peak_day)
  
df_daily_inf_filtered_sims_term2 <- ds_daily_inf_term2 %>%
  filter(
    scenario %in% c("Central_cover60_LAIV_5-18yo"),
  ) %>%
  select(c(scenario, simulation_index, run_nr, age_group, horizon, value)) %>%
  group_by(scenario, simulation_index, run_nr, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect() %>%
  group_by(scenario, simulation_index, run_nr) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    peak_day = horizon
    ) %>%
  mutate(type = "term2") %>%
  select(type, peak_day)

df_daily_inf_filtered_sims <- rbind(df_daily_inf_filtered_sims_term1, df_daily_inf_filtered_sims_term2)

df_daily_inf_filtered_summary <- df_daily_inf_filtered_sims %>%
  group_by(type) %>%
  summarise(
    median_tot = median(peak_day),
    upr95_tot = quantile(peak_day, 0.975),
    lwr95_tot = quantile(peak_day, 0.025),
    upr50_tot = quantile(peak_day, 0.75),
    lwr50_tot = quantile(peak_day, 0.25),
    .groups = "drop"
  )

fig_peak_per_term_vacc <- ggplot(
  df_daily_inf_filtered_sims,
  aes(
    x = type,
    y = peak_day,
    color = type
  )
) +
  
  # Individual model particles
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.35,
      jitter.height = 0,
      dodge.width = 0.75
    ),
    alpha = 0.2,
    size = 1.0,
    shape = 16
  ) +
  
  geom_pointrange(
    data = df_daily_inf_filtered_summary ,
    aes(
      x = type,
      y = median_tot,
      ymin = lwr50_tot,
      ymax = upr50_tot,
      group = type
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.65,
    fatten = 2
  ) +
  
  scale_color_brewer(
    "Vaccination term",
    palette = "Dark2",
    # labels = c(
    #   "AIHW" = "Baseline (AIHW)",
    #   "FluCAN - high" = "High (FluCAN adjusted)",
    #   "FluCAN - low" = "Low (FluCAN adjusted)"
    # ),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, size = 2)
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult=c(0.02, 0.01)),
    labels = scales::label_comma()
  ) +
  
  labs(
    x = "Vaccination term",
    y = "Peak day of infection incidence"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    )
  ) 

ggsave(
  filename = file.path(
    output_path,
    "peak_day_per_vacc_term.png"
  ),
  plot = fig_peak_per_term_vacc,
  width = 7.33,
  height = 5.0,
  units = "in",
  device = "png",
  dpi = 300
)

################################################################################
# END EFFECT OF VACCINATION TIMING ON EPIDEMIC
################################################################################

# ################################################################################
# FIGURE AGE FACETTED CORRELATION PLOT OF PEAK DAYS OF INFECTION INCIDENCE
# ################################################################################ 
 
source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning")
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_baseline_only <- ds_daily_inf %>%
  filter(scenario == "Status_quo") %>%
  select(simulation_index, run_nr, age_group, horizon, value) %>%
  collect() %>%
  group_by(simulation_index, run_nr, age_group) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-value)

df_baseline_only$age_group <- factor(df_baseline_only$age_group, levels = age_levels_no_all, ordered = TRUE)

df_pairs <- df_baseline_only %>%
  select(
    simulation_index,
    run_nr,
    age_x = age_group,
    peak_x = horizon
  ) %>%
  inner_join(
    df_baseline_only %>%
      select(
        simulation_index,
        run_nr,
        age_y = age_group,
        peak_y = horizon
      ),
    by = c("simulation_index", "run_nr"),
    relationship = "many-to-many"
  ) %>%
  filter(age_y > age_x)

# fig_peak_pairs <- ggplot(
#   df_pairs,
#   aes(
#     x = peak_x,
#     y = peak_y
#   )
# ) +
#   
#   geom_abline(
#     data = df_pairs %>% distinct(age_x, age_y),
#     aes(intercept = 0, slope = 1),
#     linewidth = 0.3,
#     linetype = "dashed",
#     colour = "black",
#     inherit.aes = FALSE
#   ) +
#   
#   geom_point(
#     alpha = 0.2,
#     size = 0.7,
#     shape = 16,
#     colour = "royalblue"
#   ) +
#   
#   facet_grid(
#     rows = vars(age_y),
#     cols = vars(age_x),
#     # drop = FALSE
#   ) +
#   
#   scale_x_continuous(
#     breaks = scales::breaks_width(100)
#   ) +
#   
#   labs(
#     x = "Peak day of infections",
#     y = "Peak day of infections"
#   ) +
#   
#   common_theme +
#   
#   theme(
#     # Show both horizontal and vertical major gridlines
#     panel.grid.major.x = element_line(
#       colour = "grey85",
#       linewidth = 0.3
#     ),
#     panel.grid.major.y = element_line(
#       colour = "grey85",
#       linewidth = 0.3
#     ),
#     
#     # Smaller axis tick labels
#     axis.text = element_text(size = 8)
#   )
# 
# ggsave(
#   filename = file.path(
#     output_path,
#     "fig_peak_pairs_age_strat.png"
#   ),
#   plot = fig_peak_pairs,
#   width = 7.33,
#   height = 8.9,
#   units = "in",
#   device = "png",
#   dpi = 300
# )
# 
# 

age_levels <- age_levels_no_all

n_age <- length(age_levels)
plot_list <- list()

for (row in 2:n_age) {
  
  row_plots <- list()
  
  for (col in 1:(n_age - 1)) {
    
    if (col < row) {
      
      age_x_this <- age_levels[col]
      age_y_this <- age_levels[row]
      
      df_this <- df_pairs %>%
        filter(
          age_x == age_x_this,
          age_y == age_y_this
        )
      
      p <- ggplot(
        df_this,
        aes(
          x = peak_x,
          y = peak_y
        )
      ) +
        
        geom_abline(
          intercept = 0,
          slope = 1,
          linewidth = 0.3,
          linetype = "dashed"
        ) +
        
        geom_point(
          alpha = 0.2,
          size = 0.6,
          shape = 16,
          colour = "royalblue"
        ) +
        
        scale_x_continuous(
          breaks = seq(0, 400, by = 100)
        ) +
        
        common_theme +
        
        theme(
          panel.grid.major.x = element_line(
            colour = "grey85",
            linewidth = 0.3
          ),
          panel.grid.major.y = element_line(
            colour = "grey85",
            linewidth = 0.3
          ),
          axis.text = element_text(size = 7),
          
          # Suppress axis titles on individual panels
          axis.title = element_blank()
        )
      
      # Only show x-axis tick labels on bottom row
      if (row != n_age) {
        p <- p +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
      }
      
      # Only show y-axis tick labels on first column
      if (col != 1) {
        p <- p +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
      }
      
      # Age-group label along bottom
      if (row == n_age) {
        p <- p +
          labs(x = age_x_this) +
          theme(
            axis.title.x = element_text(size = 8)
          )
      }
      
      # Age-group label along left
      if (col == 1) {
        p <- p +
          labs(y = age_y_this) +
          theme(
            axis.title.y = element_text(size = 8)
          )
      }
      
    } else {
      
      # Completely blank space above diagonal
      p <- plot_spacer()
      
    }
    
    row_plots[[col]] <- p
    
  }
  
  plot_list[[row - 1]] <- wrap_plots(
    row_plots,
    nrow = 1
  )
}

fig_peak_pairs <- wrap_plots(
  plot_list,
  ncol = 1
) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(5, 5, 5, 5)
    )
  )

ggsave(
  filename = file.path(
    output_path,
    "fig_peak_pairs_age_strat_v2.png"
  ),
  plot = fig_peak_pairs,
  width = 7.33,
  height = 8.9,
  units = "in",
  device = "png",
  dpi = 300
)

# ################################################################################
# END FIGURE AGE FACETTED CORRELATION PLOT OF PEAK DAYS OF INFECTION INCIDENCE
# ################################################################################ 

# ################################################################################
# FIGURE VACCINE ROLLOUT TIMING
# ################################################################################ 


iiv_rates <- c(
  0.00266946, 0.00334557, 0.00737364, 0.03162889, 0.07259051, 0.1482384,
  0.24145878, 0.36451686, 0.47168377, 0.56855491, 0.65137032, 0.75339126,
  0.81009063, 0.84931351, 0.88319599, 0.91812563, 0.93545698, 0.95122565,
  0.96457128, 0.97785923, 0.9846138, 0.99014375, 0.99475859, 1.0
)

laiv_rates <- c(
  0.16666667, 0.33333333, 0.5,
  0.66666667, 0.83333333, 1
)

week_end_dates <- as.Date(c(
  "2026-03-07", "2026-03-14", "2026-03-21", "2026-03-31",
  "2026-04-07", "2026-04-14", "2026-04-21", "2026-04-30",
  "2026-05-07", "2026-05-14", "2026-05-21", "2026-05-31",
  "2026-06-07", "2026-06-14", "2026-06-21", "2026-06-30",
  "2026-07-07", "2026-07-14", "2026-07-21", "2026-07-31",
  "2026-08-07", "2026-08-14", "2026-08-21", "2026-08-31"
))


df_vaccination <- bind_rows(
  data.frame(
    week_end = week_end_dates,
    cumulative = iiv_rates,
    vaccine = "IIV"
  ),
  
  data.frame(
    week_end = week_end_dates,
    cumulative = c(
      laiv_rates,
      rep(1, length(iiv_rates) - length(laiv_rates))
    ),
    vaccine = "LAIV"
  )
)

fig_vaccination_rollout <- ggplot(df_vaccination, aes(x = week_end, y = cumulative, colour = vaccine)) +
  geom_point(alpha = 0.7, size = 1.0, shape = 16) +
  geom_line() +
  scale_color_brewer("Vaccination type",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(alpha = 1, size = 2)
                     ))+
  scale_x_date(
        breaks = week_end_dates,
        date_labels = "%b %d"
      ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult=c(0.02, 0.05)),
    breaks= seq(0, 1, by = 0.2),
    labels = scales::percent
  )+
  labs(
    x = "Week ending on",
    y = "Campaign progression (%)"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
  ) 

ggsave(filename = file.path(output_path,"fig_vaccination_rollout.png"), plot = fig_vaccination_rollout, height = 4.5, width = 7.33, units = "in", device = "png", dpi = 300)


# ################################################################################
# END FIGURE VACCINE ROLLOUT TIMING
# ################################################################################ 
# 
# ################################################################################
# FIGURE TIMING PEAK SCHOOL-AGE CHILDREN VS REST
# ################################################################################ 


source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning")
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_timing_peaks <- ds_daily_inf %>%
  filter(
    scenario == "Status_quo"
  ) %>%
  select(simulation_index, run_nr, age_group, horizon, value) %>%
  collect()

df_filtered_timing_peaks_school_age <- df_timing_peaks %>%
  filter(age_group %in% c("5-11", "12-17")) %>%
  group_by(simulation_index, run_nr, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(simulation_index, run_nr) %>%
  slice_max(
    order_by = value,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(age_category = "school age") %>%
  select(age_category, peak_day = horizon)

df_filtered_timing_peaks_non_school_age <- df_timing_peaks %>%
  filter(!age_group %in% c("5-11", "12-17")) %>%
  group_by(simulation_index, run_nr, horizon) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(simulation_index, run_nr) %>%
  slice_max(
    order_by = value,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(age_category = "non school age") %>%
  select(age_category, peak_day = horizon)

df_filtered_timing_peaks <- rbind(df_filtered_timing_peaks_school_age, df_filtered_timing_peaks_non_school_age)

df_filtered_timing_peaks$age_category <- factor(df_filtered_timing_peaks$age_category, levels = c("school age", "non school age"))

df_filtered_timing_peaks_summary <- df_filtered_timing_peaks %>%
  group_by(age_category) %>%
  summarise(
    median_tot = median(peak_day),
    upr95_per = quantile(peak_day, 0.975),
    lwr95_per = quantile(peak_day, 0.025),
    upr50_per = quantile(peak_day, 0.75),
    lwr50_per = quantile(peak_day, 0.25),
    .groups = "drop"
  )

median_breaks <- df_filtered_timing_peaks_summary$median_tot

fig_timing_peaks <- ggplot(df_filtered_timing_peaks, aes(x=age_category, y= peak_day, color=age_category))+
  geom_point(position = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.75), alpha=0.2, size = 1.0, shape=16)+
  geom_pointrange(
    data = df_filtered_timing_peaks_summary,
    aes(x=age_category, y=median_tot, ymin=lwr50_per, ymax=upr50_per, group=age_category),
    position = position_dodge(width = 0.75), inherit.aes = FALSE, colour = "black", linewidth = 0.65 , fatten = 2)+
  scale_color_brewer("School age status",palette="Dark2",
                     guide = guide_legend(
                       override.aes = list(alpha = 1, size = 2)
                     ))+
  scale_y_continuous(
    limits = c(0, 365),
    expand = expansion(mult=c(0.02, 0.01)),
    breaks = sort(unique(c(
      c(0, 100, 200, 300)
    ))),
    labels = scales::label_number(accuracy = 1)
  )+
  labs(
    x = "School age status",
    y = "Day of peak infection incidence"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "black"
  ) +
  
  geom_segment(
    x = 0,
    xend = 1,
    y = df_filtered_timing_peaks_summary$median_tot[1],
    yend = df_filtered_timing_peaks_summary$median_tot[1],
    linetype = "dashed",
    linewidth = 0.35,
    colour = "black"
  ) +
  
  geom_segment(
    x = 0,
    xend = 2,
    y = df_filtered_timing_peaks_summary$median_tot[2],
    yend = df_filtered_timing_peaks_summary$median_tot[2],
    linetype = "dashed",
    linewidth = 0.35,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = 0.55,
    y = median_breaks[1] - 11,
    label = round(median_breaks[1]),
    hjust = 1
  ) +
  
  annotate(
    "text",
    x = 0.55,
    y = median_breaks[2] + 11,
    label = round(median_breaks[2]),
    hjust = 1
  ) +
  
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    )
  ) 

ggsave(filename = file.path(output_path,"fig_timing_peaks.png"), plot = fig_timing_peaks, height = 4.5, width = 7.33, units = "in", device = "png", dpi = 300)

# ################################################################################
# END FIGURE TIMING PEAK SCHOOL-AGE CHILDREN VS REST
# ################################################################################ 

# ################################################################################
# FIGURE CONTACT MATRIX
# ################################################################################ 

age_groups <- c("<1", "1–2", "3–4", "5–11",
                "12–17", "18–64", "65–79", "80+")

contact_matrix <-t(matrix(
  c(  
    0.375779488,	0.221425972,	0.117743201,	0.052502325,	0.02248264,	0.093521173,	0.007517928,	0.001580057,
    0.224587717,	0.68020258,	0.275817706,	0.091685117,	0.032435929,	0.103449184,	0.015630453,	0.00444946,
    0.368313203,	0.850640683,	2.837925261,	0.707209083,	0.14866367,	0.325028629,	0.117721982,	0.048703109,
    0.405668348,	0.698447857,	1.746863481,	8.467514786,	1.801643719,	0.677721934,	0.347731176,	0.189136569,
    0.154319918,	0.219504573,	0.326210601,	1.600482745,	9.463746013,	0.878385651,	0.189544108,	0.126039369,
    5.355102935,	5.840197141,	5.949740779,	5.022464015,	7.327712821,	10.38538854,	3.609060442,	2.608662633,
    0.091321412,	0.187192551,	0.457141402,	0.546670304,	0.3354363,	0.765615499,	1.753038513,	0.728460955,
    0.006196919,	0.017204937,	0.061063042,	0.096003292,	0.072016955,	0.178674775,	0.235198769,	0.129056892
    ),
  nrow = 8,
  byrow = TRUE
))

# Convert matrix to long format
df_contact <- as.data.frame(contact_matrix) %>%
  mutate(age_from = age_groups) %>%
  pivot_longer(
    cols = -age_from,
    names_to = "age_to_index",
    values_to = "contact"
  ) %>%
  mutate(
    age_to = rep(age_groups, times = length(age_groups)),
    age_from = factor(age_from, levels = age_groups),
    age_to = factor(age_to, levels = age_groups)
  )

age_labels <- c(
  "<1"    = "<1",
  "1–2"   = "1–<3",
  "3–4"   = "3–<5",
  "5–11"  = "5–<12",
  "12–17" = "12-<18",
  "18–64" = "18-<65",
  "65–79" = "65-<80",
  "80+"   = "80+"
)

# Plot
fig_contact_matrix <- ggplot(
  df_contact,
  aes(
    x = age_to,
    y = age_from,
    fill = contact
  )
) +
  geom_tile() +
  scale_x_discrete(labels = age_labels) +
  scale_y_discrete(labels = age_labels) +
  # scale_fill_gradient(
    # low = "aliceblue",
    # high = "darkblue"
  # ) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    transform = 'sqrt'
    ) +
  labs(
    x = "Age group to",
    y = "Age group from"
  ) +
  coord_equal() +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

ggsave(filename = file.path(output_path,"Conmat_matrix.png"), plot = fig_contact_matrix, height = 4.5, width = 7.33, units = "in", device = "png", dpi = 300)

# ################################################################################
# END FIGURE CONTACT MATRIX
# ################################################################################ 

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

source_file_daily_inf <- file.path(ROOT, paste("R outputs", source_file_type, vaccination_term), "dat_uom_infection_with_waning")
ds_daily_inf <- open_dataset(source_file_daily_inf)

df_temp <- ds_daily_inf %>%
  filter(
    scenario == "Status_quo"
  ) %>%
  group_by(simulation_index, run_nr, age_group) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  filter(age_group %in% c("5-11", "12-17"))%>%
  collect() %>%
  left_join(ages, by = "age_group") %>%
  mutate(ipht = 100000 * value / pop_size)

df_temp_5_11 <- df_temp %>%
  filter(age_group == "5-11")

min(df_temp_5_11$ipht)

max(df_temp_5_11$ipht)

median(df_temp_5_11$ipht)

df_temp_12_17 <- df_temp %>%
  filter(age_group == "12-17")

min(df_temp_12_17$ipht)

max(df_temp_12_17$ipht)

median(df_temp_12_17$ipht)


####################################################################################


# df_R0 <- read_csv(file.path(ROOT, "ABC outputs with waning", "list_R0.csv")) %>%
#   select(R0)
# 
# min(df_R0$R0)
# 
# max(df_R0$R0)
# 
# median(df_R0$R0)



















