
source_file_type <- "with waning"
source_file_formatted <- gsub("\\s+", "_", source_file_type)

if (.Platform$OS.type == "windows"){
  setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Code")
  output_root <- file.path(getwd(), "Scenario-hub-visualisations", "data")
  source_file_dir <- file.path(paste("Forward projection",source_file_type), paste("df_incid_", source_file_formatted, ".parquet", sep = ""))
} else{
  setwd("/home/ubuntu/R_code")
  output_root <- "/pvol"
  source_file_path <- file.path(output_root, paste("Forward projection", source_file_type), paste("df_incid_", source_file_formatted, ".parquet", sep="") )
}

library(tidyverse)
library(arrow)
library(tictoc)

options(readr.show_col_types = FALSE)

tic()
tic()

clean_age_group <- function(age_group){
  dplyr::case_when(
    age_group %in% c("0-1")    ~ "<1",
    age_group %in% c("1-<5")   ~ "1-4",
    age_group %in% c("5-<12")  ~ "5-11",
    age_group %in% c("12-<18") ~ "12-17",
    age_group %in% c("18-<65") ~ "18-64",
    age_group %in% c("65-<80") ~ "65-79",
    age_group %in% c(">80")    ~ "80+",
    TRUE ~ age_group)
}

age_to_index <- function(age_group){
  case_when(
    age_group == "<1" ~ 1,
    age_group == "1" ~ 2,
    age_group == "2-4" ~ 2,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )
}

add_age_and_scenario_labels <- function(df){
  df %>%
    mutate(
      age_group = case_when(
        age_index == 1 ~ "<1",
        age_index == 2 ~ "1-4",
        age_index == 3 ~ "5-11",
        age_index == 4 ~ "12-17",
        age_index == 5 ~ "18-64",
        age_index == 6 ~ "65-79",
        age_index == 7 ~ "80+"
      ),
      scenario = case_when(
        scenario_index == 1 ~ "Status_quo",
        scenario_index == 2 ~ "Pessimistic_LAIV_5-12yo",
        scenario_index == 3 ~ "Pessimistic_LAIV_5-18yo",
        scenario_index == 4 ~ "Mid_LAIV_5-12yo",
        scenario_index == 5 ~ "Mid_LAIV_5-18yo",
        scenario_index == 6 ~ "Optimistic_LAIV_5-12yo",
        scenario_index == 7 ~ "Optimistic_LAIV_5-18yo",
        scenario_index == 8 ~ "Mid_LAIV_2-5yo",
        scenario_index == 9 ~ "Mid_LAIV_2-12yo",
        scenario_index == 10 ~ "Mid_LAIV_2-18yo",
        scenario_index == 11 ~ "Mid_Low_LAIV_5-12yo",
        scenario_index == 12 ~ "Mid_Low_LAIV_5-18yo",
        scenario_index == 13 ~ "Mid_Low_LAIV_2-12yo",
        scenario_index == 14 ~ "Mid_Low_LAIV_2-18yo",
        scenario_index == 15 ~ "Mid_Extra_Low_LAIV_5-12yo",
        scenario_index == 16 ~ "Mid_Extra_Low_LAIV_5-18yo"
      )
    ) %>%
    select(-age_index, -scenario_index)
}

# hospitalisations data
df_hosp <- read_csv(file.path("data", "all_aihw_hosps_raw_2023.csv")) %>%
  select(-Population, -Separations) %>%
  filter(Age != "All") %>%
  rename(age_group = Age, hpht = Rate) %>% 
  mutate(
    age_group = clean_age_group(age_group),
    age_index = age_to_index(age_group)) %>%
  select(-age_group)

# cases data
df_case <- read_csv(file.path("data", "cases_2023_by_100000_age.csv")) %>%
  select(-1) %>%
  rename(cpht = cases_per_100000_age) %>%
  mutate(
    age_group = clean_age_group(age_group),
    age_index = age_to_index(age_group )) %>%
  select(-age_group)

# fatalities data
df_fatal <- read_csv(file.path("data", "influenza_deaths_cleaned.csv")) %>%
  rename(fpht = Rate_per_100K, age_group = Age) %>%
  select(-c(population, Deaths)) %>%
  filter(age_group != "All") %>%
  mutate(
    age_group = clean_age_group(age_group),
    age_index = age_to_index(age_group)
  ) %>%
  select(-age_group)

# Load age group sizes from simulation to rescale to 100k per age group
ages <- read_csv(file.path("data", "age_pops.csv")) %>%
  mutate(`1-4` = `1` + `2-4`) %>%
  select(-`1`, -`2-4`) %>%
  pivot_longer(cols = everything(), names_to = "age_group",
               values_to = "pop_size")%>% 
  mutate(age_index = age_to_index(age_group)) %>%
  select(age_index, pop_size)

ds_case  <- Table$create(df_case)
ds_hosp  <- Table$create(df_hosp)
ds_fatal <- Table$create(df_fatal)
ds_ages  <- Table$create(ages)

# replace the string valued columns by integers to decrease memory usage of further manipulations

ds_full <- open_dataset(source_file_path) %>%
  mutate(
    scenario_index = case_when(
      scenario == "Status quo" ~ 1,
      scenario == "Pessimistic_LAIV_5-12yo" ~ 2,
      scenario == "Pessimistic_LAIV_5-18yo" ~ 3,
      scenario == "Mid_LAIV_5-12yo" ~ 4,
      scenario == "Mid_LAIV_5-18yo" ~ 5,
      scenario == "Optimistic_LAIV_5-12yo" ~ 6,
      scenario == "Optimistic_LAIV_5-18yo" ~ 7,
      scenario == "Mid_LAIV_2-5yo" ~ 8,
      scenario == "Mid_LAIV_2-12yo" ~ 9,
      scenario == "Mid_LAIV_2-18yo" ~ 10,
      scenario == "Mid_Low_LAIV_5-12yo" ~ 11,
      scenario == "Mid_Low_LAIV_5-18yo" ~ 12,
      scenario == "Mid_Low_LAIV_2-12yo" ~ 13,
      scenario == "Mid_Low_LAIV_2-18yo" ~ 14,
      scenario == "Mid_Extra_Low_LAIV_5-12yo" ~ 15,
      scenario == "Mid_Extra_Low_LAIV_5-18yo" ~ 16
      ),
    target_index = case_when(
      target == "infection_incidence_vacc" ~ 1,
      target == "infection_incidence_unvacc" ~ 2
    ),
    # give both age groups "1" and "2-4" the code 2 in order to aggregate them below
    age_index = case_when(
      age_group == "<1" ~ 1,
      age_group == "1" ~ 2,
      age_group == "2-4" ~ 2,
      age_group == "5-11" ~ 3,
      age_group == "12-17" ~ 4,
      age_group == "18-64" ~ 5,
      age_group == "65-79" ~ 6,
      age_group == "80+" ~ 7
    )) %>% 
  select(scenario_index, simulation_index, age_index, horizon, run_nr, target_index, value)

# dataset with age groups 1 and 2-4 aggregated into 1-4
tmp_collapsed_path <- file.path(output_root, "tmp_full_collapsed")

if (dir.exists(tmp_collapsed_path)) {
  unlink(tmp_collapsed_path, recursive = TRUE)
}

scenario_vals <- ds_full %>%
  select(scenario_index) %>%
  distinct() %>%
  arrange(scenario_index) %>%
  collect() %>%
  pull(scenario_index)

nr_scenarios <- length(scenario_vals)

for (s in scenario_vals) {
  
  message("Processing scenario ", s, " out of ", nr_scenarios, ", collapsing relevant age groups")
  
  ds_full %>%
    filter(scenario_index == s) %>%
    group_by(
      scenario_index, simulation_index, age_index,
      horizon, run_nr, target_index
    ) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    write_dataset(
      tmp_collapsed_path,
      format = "parquet",
      partitioning = "scenario_index",
      existing_data_behavior = "delete_matching"
    )
  invisible(gc())
}

ds_full_collapsed <- open_dataset(tmp_collapsed_path)

car_ihr <- ds_full_collapsed %>%
  filter(scenario_index == 1L) %>%
  group_by(scenario_index, simulation_index, age_index, run_nr, target_index) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(scenario_index, simulation_index, age_index, target_index) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  left_join(ds_ages, by = "age_index") %>%
  mutate(ipht = 1e5 * value / pop_size) %>%
  select(-value, -pop_size) %>%
  collect()

car_ihr_vacc <- car_ihr %>%
  filter(target_index == 1L) %>%
  rename(inf_inc_vacc = ipht) %>%
  select(-target_index)

car_ihr <- car_ihr %>%
  filter(target_index == 2L) %>%
  rename(inf_inc_unvacc = ipht) %>%
  select(-target_index) %>%
  left_join(
    car_ihr_vacc,
    by = c("scenario_index", "simulation_index", "age_index")
  ) %>%
  left_join(df_case, by = "age_index") %>%
  left_join(df_hosp, by = "age_index") %>%
  left_join(df_fatal, by = "age_index") %>%
  mutate(
    car = cpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    ihr = hpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    ifr = fpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc)
  ) %>%
  select(simulation_index, age_index, car, ihr, ifr)

car_ihr_labelled <- car_ihr %>%
  mutate(
    age_group = case_when(
      age_index == 1 ~ "<1",
      age_index == 2 ~ "1-4",
      age_index == 3 ~ "5-11",
      age_index == 4 ~ "12-17",
      age_index == 5 ~ "18-64",
      age_index == 6 ~ "65-79",
      age_index == 7 ~ "80+"
    )
  ) %>%
  select(-age_index)

output_path <- file.path(output_root, paste("R outputs", source_file_type) )

if (dir.exists(output_path)) {
  unlink(output_path, recursive = TRUE)
} 

dir.create(output_path, recursive = TRUE)

write_parquet(car_ihr_labelled, file.path(output_path, paste("car_ihr_ifr_", source_file_formatted, ".parquet", sep = "")))

rm(car_ihr_vacc, car_ihr_labelled)
invisible(gc())

time_elapsed <- toc(log = TRUE, quiet = TRUE)
cat("Creation of car-ihr file took", round(time_elapsed$toc - time_elapsed$tic, 2), "seconds.")

################################################################

tmp_wide_path <- file.path("data", "tmp_full_wide")

if (dir.exists(tmp_wide_path)) {
  unlink(tmp_wide_path, recursive = TRUE)
}

for (s in scenario_vals) {
  
  message("Processing scenario ", s, " out of ", nr_scenarios, ", widening vacc - unvacc columns")
  
  ds_full_collapsed %>%
    filter(scenario_index == s) %>%
    group_by(
      scenario_index, simulation_index, age_index, horizon, run_nr
    ) %>%
    summarise(
      inf_inc_unvacc = sum(if_else(target_index == 2L, value, 0)),
      inf_inc_vacc   = sum(if_else(target_index == 1L, value, 0)),
      .groups = "drop"
    ) %>%
    write_dataset(
      tmp_wide_path,
      format = "parquet",
      partitioning = "scenario_index",
      existing_data_behavior = "delete_matching"
    )
  
  invisible(gc())
}

ds_full_wide <- open_dataset(tmp_wide_path)

ds_car_ihr <- Table$create(
  car_ihr %>%
    mutate(
      simulation_index = as.integer(simulation_index),
      age_index = as.integer(age_index)
    )
)

ds_full_enriched <- ds_full_wide %>%
  mutate(
    simulation_index = as.integer(simulation_index),
    age_index = as.integer(age_index)
  ) %>%
  left_join(ds_car_ihr, by = c("simulation_index", "age_index"))%>%
  mutate(
    infection_incidence = inf_inc_unvacc + inf_inc_vacc,
    disease_incidence = car * (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    admission_incidence = ihr * (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    fatality_incidence = ifr * (inf_inc_unvacc + 0.8 * inf_inc_vacc)
  ) %>%
  select(
    scenario_index, simulation_index, age_index, horizon, run_nr,
    infection_incidence, disease_incidence,
    admission_incidence, fatality_incidence
  )

write_target <- function(ds, value_col, target_name, out_path) {
  ds %>%
    transmute(
      scenario_index,
      simulation_index,
      age_index,
      horizon,
      run_nr,
      value = !!sym(value_col),
      target = target_name,
      setting = "temperate",
      team = "UoM"
    ) %>%
    add_age_and_scenario_labels() %>%
    write_parquet(out_path)
  
  invisible(gc())
}

write_target(
  ds_full_enriched,
  "infection_incidence",
  "infection_incidence",
  file.path(output_path, paste("dat_uom_infection_", source_file_formatted, ".parquet", sep = ""))
)

write_target(
  ds_full_enriched,
  "disease_incidence",
  "disease_incidence",
  file.path(output_path, paste("dat_uom_disease_", source_file_formatted, ".parquet", sep = ""))
)

write_target(
  ds_full_enriched,
  "admission_incidence",
  "admission_incidence",
  file.path(output_path, paste("dat_uom_admission_", source_file_formatted, ".parquet", sep = ""))
)

write_target(
  ds_full_enriched,
  "fatality_incidence",
  "fatality_incidence",
  file.path(output_path, paste("dat_uom_fatality_", source_file_formatted, ".parquet", sep = ""))
)

unlink(tmp_collapsed_path, recursive = TRUE)
unlink(tmp_wide_path, recursive = TRUE)

time_elapsed <- toc(log = TRUE, quiet = TRUE)
cat("Full run took", round(time_elapsed$toc - time_elapsed$tic, 2), "seconds.\n")

# Just some testing stuff

# names(ds_full)

# ds_full %>%
#   summarise(n = n()) %>%
#   collect()
# # 
# ds_full %>%
#   distinct(simulation_index) %>%
#   collect()
# 
# ds_full %>%
#   summarise(n = n_distinct(simulation_index)) %>%
#   collect()

# ds_full %>%
#   mutate(scenario = as.character(scenario)) %>%
#   distinct(scenario) %>%
#   collect()

