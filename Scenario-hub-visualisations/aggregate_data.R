source_file_type <- "with waning"
vaccination_term <- "term2 vaccination"

library(arrow)
library(dplyr)

source_file_formatted <- gsub("\\s+", "_", source_file_type)
# target_list <- c("infection")
# target_list <- c("admission_worst_flucan")
target_list <- c("infection", "disease","admission", "admission_best_flucan", "admission_worst_flucan")

for (target in target_list){
  
  cat("Aggregating", target, "results.\n")
  
  if (.Platform$OS.type == "windows"){
    setwd("C:/Users/ebisa/Documents/University of Melbourne/Scenario modelling exercise 2025/Code")
    output_root <- file.path(getwd(), "Scenario-hub-visualisations", "data")
    source_file_path <- file.path(paste("Forward projection",source_file_type), paste("df_incid_", source_file_formatted, ".parquet", sep = ""))
  } else{
    setwd("/home/ubuntu/R_code")
    output_root <- "/pvol"
    source_file_path <- file.path(output_root, paste("R outputs", source_file_type, vaccination_term), paste("dat_uom", target ,source_file_formatted, sep="_") )
    output_path <- file.path(output_root, paste("R outputs", source_file_type, vaccination_term))
  }
  
  target_ds <- open_dataset(source_file_path)
  
  target_ar_ds <- target_ds %>%
    group_by(scenario, simulation_index, run_nr, age_group) %>%
    summarise(value = sum(value), .groups = "drop") 
  
  unlink(file.path(output_path, target_ar_ds), recursive = TRUE)
  
  write_dataset(
    target_ar_ds,
    path = file.path(output_path, paste(target ,"ar_ds", sep="_")),
    format = "parquet"
  )
  cat("Finished aggregating", target, "results.\n\n")
}
 





