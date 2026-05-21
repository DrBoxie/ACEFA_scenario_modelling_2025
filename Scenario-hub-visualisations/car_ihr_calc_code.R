
# setwd("/home/ubuntu"")
# setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Scenario-hub-visualisations")

library(tidyverse)
library(arrow)
library(tictoc)
library(DBI)
library(duckdb)

ds <- open_dataset("data/dat_uom_inf_and_dis.parquet")

names(ds)

ds %>%
  summarise(n = n()) %>%
  collect()

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

# hospitalisations data
df_hosp <- read_csv(file.path("data", "all_aihw_hosps_raw_2023.csv")) %>%
  select(-Population, -Separations) %>%
  filter(Age != "All") %>%
  rename(age_group = Age, hpht = Rate) %>%
  mutate(age_group = clean_age_group(age_group))

# cases data
df_case <- read_csv(file.path("data", "cases_2023_by_100000_age.csv")) %>%
  select(-1) %>%
  rename(cpht = cases_per_100000_age) %>%
  mutate(age_group = clean_age_group(age_group))

df_fatal <- read_csv(file.path("data", "influenza_deaths_cleaned.csv")) %>%
  rename(fpht = Rate_per_100K, age_group = Age) %>%
  select(-c(population, Deaths)) %>%
  filter(age_group != "All")

# Load age group sizes from simulation to rescale to 100k per age group
ages <- read_csv(file.path("data", "age_pops.csv")) %>%
  mutate(`1-4` = `1` + `2-4`) %>%
  select(-`1`, -`2-4`) %>%
  pivot_longer(cols = everything(), names_to = "age_group",
               values_to = "pop_size")%>% 
  mutate(age_index = case_when(
    age_group == "<1" ~ 1,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )) %>%
  select(-age_group)

# transform the dataframes above into a format compatible with lazy querying
ds_hosp <- df_hosp %>%
  mutate(age_index = case_when(
    age_group == "<1" ~ 1,
    age_group == "1" ~ 2,
    age_group == "2-4" ~ 2,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )) %>% 
  select(-age_group) %>%
  Table$create()

ds_case <- df_case %>%
  mutate(age_index = case_when(
    age_group == "<1" ~ 1,
    age_group == "1" ~ 2,
    age_group == "2-4" ~ 2,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )) %>% 
  select(-age_group) %>%
  Table$create()

ds_fatal <- df_fatal %>%
  mutate(age_index = case_when(
    age_group == "<1" ~ 1,
    age_group == "1" ~ 2,
    age_group == "2-4" ~ 2,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )) %>% 
  select(-age_group) %>%
  Table$create()
  
ds_ages <- Table$create(ages)

# replace the string valued columns by integers to decrease memory usage of further manipulations
ds_full <- open_dataset(file.path("data", "df_incid.parquet")) %>%
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
      scenario == "Mid_LAIV_2-18yo" ~ 10
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
  select(-c(age_group, target, scenario, setting))

# dataset with age groups 1 and 2-4 aggregated into 1-4
ds_full_collapsed <- union_all(
  # Collapse age groups 1 and 2-4 into 1-4
  ds_full %>% filter(age_index == 2) %>%
    group_by(scenario_index, simulation_index, age_index, horizon, run_nr, target_index) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    
    # Reorder to match original ds_full column ordering
    select(simulation_index, horizon, run_nr, value, scenario_index, target_index, age_index),
  
  # Collect all other age groups
  ds_full %>% filter(age_index != 2)
)

ds_baseline <- ds_full_collapsed %>%
  filter(scenario_index == 1) %>%                                                       # calculation of car, ihr, and ifr, are all done on the baseline scenario
  group_by(scenario_index, simulation_index, age_index, run_nr, target_index) %>%       # group and then sum across all days to get total AR for each group
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(scenario_index, simulation_index, age_index, target_index) %>%               # group and then take the mean across all runs per scenario-particle combo
  summarise(value = mean(value), .groups = "drop") %>%
  left_join(ds_ages, by = "age_index") %>%                                              # add population size needed for per 100k calculation 
  mutate(ipht = 100000 * value / pop_size) %>%
  select(-value, -pop_size)
  
ds_baseline_vacc <- ds_baseline %>%
  filter(target_index == 1) %>%
  rename(inf_inc_vacc = ipht) %>%
  select(-target_index)

ds_baseline <- ds_baseline %>%
  filter(target_index == 2) %>%
  rename(inf_inc_unvacc = ipht) %>%
  select(-target_index) %>%
  left_join(ds_baseline_vacc, by = c("scenario_index", "simulation_index", "age_index")) %>%
  left_join(ds_case, by = "age_index") %>%
  left_join(ds_hosp, by = "age_index") %>%
  left_join(ds_fatal, by = "age_index") %>%
  mutate(
    car = cpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    ihr = hpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    ifr = fpht / (inf_inc_unvacc + 0.8 * inf_inc_vacc),
  )

car_ihr <- ds_baseline %>%
  select(simulation_index, age_index, car, ihr, ifr) %>%
  mutate(age_group = case_when(
    age_index == 1 ~ "<1",
    age_index == 2 ~ "1-4",
    age_index == 3 ~ "5-11",
    age_index == 4 ~ "12-17",
    age_index == 5 ~ "18-64",
    age_index == 6 ~ "65-79",
    age_index == 7 ~ "80+",
  )) %>%
  select(-age_index)

write_parquet(car_ihr, file.path("data", "car_ihr.parquet"))

time_elapsed <- toc(log = TRUE, quiet = TRUE)
cat("Creation of car-ihr file took", round(time_elapsed$toc - time_elapsed$tic, 2), "seconds.")

################################################################

car_ihr <- car_ihr %>%
  mutate( age_index = case_when(
    age_group == "<1" ~ 1,
    age_group == "1-4" ~ 2,
    age_group == "5-11" ~ 3,
    age_group == "12-17" ~ 4,
    age_group == "18-64" ~ 5,
    age_group == "65-79" ~ 6,
    age_group == "80+" ~ 7
  )) %>% 
  select(-age_group)

ds_full_wide <- ds_full_collapsed %>%
  group_by(scenario_index, simulation_index, age_index, horizon, run_nr) %>%
  summarise(
    inf_inc_unvacc = sum(if_else(target_index == 2, value, 0)),
    inf_inc_vacc   = sum(if_else(target_index == 1, value, 0)),
    .groups = "drop"
  )

ds_full_enriched <- ds_full_wide %>%
  left_join(car_ihr, by = c("simulation_index", "age_index")) %>%
  mutate(
    infection_incidence = inf_inc_unvacc + inf_inc_vacc,
    disease_incidence = car * (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    admission_incidence = ihr * (inf_inc_unvacc + 0.8 * inf_inc_vacc),
    fatality_incidence = ifr * (inf_inc_unvacc + 0.8 * inf_inc_vacc)
  ) %>%
  select(-c(inf_inc_unvacc, inf_inc_vacc, car, ihr, ifr))

ds_full_final_inf_and_dis <- union_all(
  ds_full_enriched %>%
    transmute(scenario_index, simulation_index, age_index, horizon, run_nr, value = infection_incidence, target = "infection_incidence", setting = "temperate", team = "UoM"),
  ds_full_enriched %>%
    transmute(scenario_index, simulation_index, age_index, horizon, run_nr, value = disease_incidence, target = "disease_incidence", setting = "temperate", team = "UoM")
)

ds_full_final_hosp_and_fatal <- union_all(
  ds_full_enriched %>%
    transmute(scenario_index, simulation_index, age_index, horizon, run_nr, value = admission_incidence, target = "admission_incidence", setting = "temperate", team = "UoM"),
  ds_full_enriched %>%
    transmute(scenario_index, simulation_index, age_index, horizon, run_nr, value = fatality_incidence, target = "fatality_incidence", setting = "temperate", team = "UoM")
)

# function to replace the age_index and simulation_index values back to their actual string values
add_age_and_scenario_labels <- function(df) {
  df %>%
    mutate(age_group = case_when(
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
      scenario_index == 10 ~ "Mid_LAIV_2-18yo"
    )) %>%
    select(-c(age_index, scenario_index))  
}

ds_full_final_inf_and_dis <- add_age_and_scenario_labels(ds_full_final_inf_and_dis)
ds_full_final_hosp_and_fatal <- add_age_and_scenario_labels(ds_full_final_hosp_and_fatal)

write_parquet(ds_full_final_inf_and_dis, file.path("data", "dat_uom_inf_and_dis.parquet"))
write_parquet(ds_full_final_hosp_and_fatal, file.path("data", "dat_uom_hosp_and_fatal.parquet"))

con <- dbConnect(duckdb())
dbExecute(con, "
  COPY (
    SELECT * FROM read_parquet('data/dat_uom_inf_and_dis.parquet')
    UNION ALL
    SELECT * FROM read_parquet('data/dat_uom_hosp_and_fatal.parquet')
  )
  TO 'data/dat_uom.parquet'
  (FORMAT PARQUET)
")

dbDisconnect(con, shutdown = TRUE)

time_elapsed <- toc(log = TRUE, quiet = TRUE)
cat("Full run took", round(time_elapsed$toc - time_elapsed$tic, 2), "seconds.")




