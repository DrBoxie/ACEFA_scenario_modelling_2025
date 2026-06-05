###############################.
## Read in the uploaded data ## 
###############################.

library(arrow)

cat("Preparing and loading UoM data.\n\n")

# The following 4 blocks calculate the AR for each particle - run - age combination
cat("Loading infections.\n\n")
dat_uom_inf <- open_dataset(file.path("data", "dat_uom_infection_no_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

dat_uom_inf_with_waning <- open_dataset(file.path("data", "dat_uom_infection_with_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

cat("Loading cases.\n\n")
dat_uom_dis <- open_dataset(file.path("data", "dat_uom_disease_no_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

dat_uom_dis_with_waning <- open_dataset(file.path("data", "dat_uom_disease_with_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

cat("Loading hospitalisations.\n\n")
dat_uom_adm <- open_dataset(file.path("data", "dat_uom_admission_no_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

dat_uom_adm_with_waning <- open_dataset(file.path("data", "dat_uom_admission_with_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

cat("Loading fatalities.\n\n")
dat_uom_fatal <- open_dataset(file.path("data", "dat_uom_fatality_no_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

dat_uom_fatal_with_waning <- open_dataset(file.path("data", "dat_uom_fatality_with_waning.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

#######################

add_all_ages <- function(df){
  df_tmp <- df %>%
    group_by(scenario, simulation_index, run_nr, target, setting, team) |>
    summarise(value = sum(value), .groups = "drop") |>
    mutate(age_group = "All")
    
  df <- df %>%
    bind_rows(df_tmp) %>%
    group_by(scenario, simulation_index, age_group, target, team) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    mutate(
      age_group = factor(age_group, levels = c("<1","1-4","5-11","12-17","18-64" ,"65-79","80+", "All"))
    )
  return(df)
}

cat("Adding All ages category to infections.\n\n")
dat_uom_inf <- add_all_ages(dat_uom_inf)
dat_uom_inf_with_waning <- add_all_ages(dat_uom_inf_with_waning)

cat("Adding All ages category to cases.\n\n")
dat_uom_dis <- add_all_ages(dat_uom_dis)
dat_uom_dis_with_waning <- add_all_ages(dat_uom_dis_with_waning)

cat("Adding All ages category to hospitalisations.\n\n")
dat_uom_adm <- add_all_ages(dat_uom_adm)
dat_uom_adm_with_waning <- add_all_ages(dat_uom_adm_with_waning)

cat("Adding All ages category to fatalities.\n\n")
dat_uom_fatal <- add_all_ages(dat_uom_fatal)
dat_uom_fatal_with_waning <- add_all_ages(dat_uom_fatal_with_waning)

# rm(dat_uom_all_ages)
# invisible(gc())

dat_uom_list <- list(dat_uom_inf, dat_uom_dis, dat_uom_adm, dat_uom_fatal)
dat_uom_list_with_waning <- list(dat_uom_inf_with_waning, dat_uom_dis_with_waning, dat_uom_adm_with_waning, dat_uom_fatal_with_waning)

cat("UoM data loaded!\n\n")


