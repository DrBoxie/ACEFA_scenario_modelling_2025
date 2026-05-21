###############################.
## Read in the uploaded data ## ------------------------------------------------
###############################.
library(arrow)

cat("Preparing and loading UoM data.\n\n")

dat_uom <- open_dataset(file.path("data", "dat_uom_SINGLE_RUN_ABC_AND_FORWARD.parquet")) %>%
  group_by(scenario, simulation_index, age_group, run_nr, target, setting, team) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  collect()

# Add rows for total infections across all ages
dat_uom_all_ages <- dat_uom |>
  group_by(scenario, simulation_index, run_nr, target, setting, team) |>
  summarise(value = sum(value), .groups = "drop") |>
  mutate(age_group = "All")

dat_uom <- dat_uom %>%
  rbind(dat_uom_all_ages) %>%
  group_by(scenario, simulation_index, age_group, target, team) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(
    age_group = factor(age_group, levels = c("<1","1-4","5-11","12-17","18-64" ,"65-79","80+", "All")),
  )

rm(dat_uom_all_ages)
invisible(gc())

cat("UoM data loaded!\n\n")


