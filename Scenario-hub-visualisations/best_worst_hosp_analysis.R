
rm(list = ls())

setwd("/home/ubuntu/R_code")
output_root <- file.path("data")

###############################################################################################################
# Loading required packages

library(tidyverse)
library(RColorBrewer)
library(ggplot2)

##############################################################################################################

unique_locations <- c("ACT","NSW","NT","QLD","SA","TAS","VIC","WA")
unique_locations_full <- c("Australian Capital Territory",
                           "New South Wales",
                           "Northern Territory",
                           "Queensland",
                           "South Australia",
                           "Tasmania",
                           "Victoria",
                           "Western Australia")
unique_seasons <- c("2021", "2022", "2023", "2024", "2025")
unique_age_groups <- c("0-1", "1-4"  , "5-11",  "12-17" ,"18-64", "65-79", ">80" )

###############################################################################################################
# Load the FluCan data 

# Influenza hospitalisations total
df_hosp_flucan <- read.csv(paste("data/flucan_dataset_2021_2025_processed", ".csv", sep=""))
df_hosp_flucan[df_hosp_flucan$agegroup1=="<1",]$agegroup1 <- "0-1"

# Population size by age group
# https://dataexplorer.abs.gov.au/vis?tm=quarterly%20population&pg=0&df[ds]=ABS_ABS_TOPICS&df[id]=ERP_Q&df[ag]=ABS&df[vs]=1.0.0&hc[Frequency]=Quarterly&pd=2014-Q2%2C&dq=1.3.97%2B8599%2BA80%2BA75%2BA70%2BA65%2BA60%2BA55%2BA50%2BA45%2BA40%2BA35%2BA30%2BA25%2BA59%2BA20%2B19%2B18%2B17%2B16%2B15%2B14%2B13%2B12%2B11%2B10%2B3%2B4%2B2%2B1%2B0%2BTOT..Q&ly[cl]=TIME_PERIOD&ly[rs]=AGE&ly[rw]=REGION&to[TIME_PERIOD]=false
ages <- read.csv(paste("data/AusAges", ".csv", sep=""))

# Age categories required for our analysis
age_groups <- list() 
age_groups[[1]] <- c("0")
#age_groups[[2]] <- c("1")
#age_groups[[3]] <- c("2","3","4")
age_groups[[2]] <- c("1","2","3","4")

age_groups[[3]] <- c("5-9","10","11")
age_groups[[4]] <- c("12","13","14","15","16","17")
age_groups[[5]] <- c("18","19",
                     "20-24","25-29",
                     "30-34","35-39",
                     "40-44","45-49",
                     "50-54","55-59",
                     "60-64")
age_groups[[6]] <- c("65-69","70-74","75-79")
age_groups[[7]] <- c("80-84", "85 and over")

###############################################################################################################
# Data is being wrangled 

# df <- data.frame()
df_age <- data.frame()

for(i in unique_locations){
  print(i)
  
  for(j in unique_seasons){

    for(k in unique_age_groups){
      
      if(j!="2025"){
        tmp_ages <- ages[ages$Region==unique_locations_full[i==unique_locations] &ages$TIME_PERIOD==paste(j,"-Q3", sep="") & ages$Age%in%age_groups[k==unique_age_groups][[1]] ,]
        
      } else{
        tmp_ages <- ages[ages$Region==unique_locations_full[i==unique_locations] &ages$TIME_PERIOD==paste("2024","-Q3", sep="") & ages$Age%in%age_groups[k==unique_age_groups][[1]] ,]
        
      }
        
      temp_hosp <- df_hosp_flucan[df_hosp_flucan$year %in% j & df_hosp_flucan$state %in% i & df_hosp_flucan$agegroup1 %in%k,]
      
      row_df <- data.frame(season = j,
                           location = i,
                           age = k,
                           pop_size = sum(tmp_ages$OBS_VALUE),
                           hosp_ar = sum(temp_hosp$admissions))
        
       
      df_age <- rbind(df_age, row_df)
    }
  }
}

################################################################################

pops <- c(sum(df_age[df_age$season=="2025" & df_age$age=="0-1",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age=="1-4",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age=="5-11",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age=="12-17",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age=="18-64",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age=="65-79",]$pop_size),
          sum(df_age[df_age$season=="2025" & df_age$age==">80",]$pop_size))

################################################################################

df_age$age <- factor(df_age$age, levels=c("0-1", "1-4", "5-11", "12-17", "18-64", "65-79", ">80"))

# Remove years 2020 and 2021
df_age <- df_age[!(df_age$season %in% c("2020","2021")),]

# Calculate hospitalisations per hundred thousand (hpht)
df_age$hpht <- df_age$hosp_ar*100000/df_age$pop_size

years_without_hosps <- unique(interaction(df_age[df_age$hpht%in%0 | is.na(df_age$hpht),]$season, df_age[df_age$hpht%in%0 | is.na(df_age$hpht),]$location ))

# df_age_input = rows with non-zero hospitalisations
df_age <- df_age[!(interaction(df_age$season, df_age$location) %in% years_without_hosps),]

hosp_aggregated_flucan <- df_age %>%
  group_by(season, location) %>%
  summarise(
    total_hosp = sum(hosp_ar),
    total_pop = sum(pop_size),
    .groups = "drop"
  ) %>%
  mutate(hpht = 100000 * total_hosp / total_pop) %>%
  arrange(hpht)

worst_hosp_year_location <- hosp_aggregated_flucan %>%
  filter(hpht == max(hpht))
worst_year <- worst_hosp_year_location$season
worst_location <- worst_hosp_year_location$location

best_hosp_year_location <- hosp_aggregated_flucan %>%
  filter(hpht == min(hpht))
best_year <- best_hosp_year_location$season
best_location <- best_hosp_year_location$location

mid_hosp_year_location <- hosp_aggregated_flucan %>%
  filter(season == "2023", location=="SA")
mid_year <- mid_hosp_year_location$season
mid_location <- mid_hosp_year_location$location

# save best and worst year dataframe

best_hosp_flucan <- df_age %>%
  filter(
    season %in% c(best_year),
    location %in% c(best_location)
  ) %>%
  mutate(type = "best")

worst_hosp_flucan <- df_age %>%
  filter(
    season %in% c(worst_year),
    location %in% c(worst_location)
  ) %>%
  mutate(type = "worst")

best_worst_hosp_flucan <- rbind(best_hosp_flucan, worst_hosp_flucan)

write.csv(best_worst_hosp_flucan, file.path("data", "best_worst_hosp_flucan.csv"), row.names = FALSE)

###############################################################################################################
# Plot best and worst years

df_age_plot <- df_age %>%
  mutate(
    highlight = case_when(
      season == best_year & location == best_location ~ "Best",
      season == worst_year & location == worst_location ~ "Worst",
      season == mid_year & location == mid_location ~ "Mid",
      TRUE ~ "Other"
    )
  )

hosp_dat <- ggplot()+
  geom_point(data=df_age_plot, aes(x=age,y=hpht, col=location, group=interaction(season, location)), position = position_dodge(width=0.2), alpha = 0.45, size = 1.8)+
  geom_line(data=df_age_plot, aes(x=age,y=hpht, col=location, group=interaction(season, location) ), position = position_dodge(width=0.2), alpha = 0.45, linewidth = 0.6)+
  # Draw the best and worst curves again on top
  geom_line(
    data = df_age_plot %>%
      filter(highlight != "Other"),
    aes(
      x = age,
      y = hpht,
      colour = location,
      group = interaction(season, location),
      linetype = highlight
    ),
    position = position_dodge(width = 0.2),
    linewidth = 1.4
  ) +
  
  geom_point(
    data = df_age_plot %>%
      filter(highlight != "Other"),
    aes(
      x = age,
      y = hpht,
      colour = location,
      group = interaction(season, location),
      shape = highlight
    ),
    position = position_dodge(width = 0.2),
    size = 3.2,
    stroke = 1.2
  ) +
  scale_y_log10()+
  scale_color_brewer("Jurisdiction",palette = "Dark2")+
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Best" = "solid",
      "Worst" = "dashed",
      "Mid" = "dotted"
    )
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "Best" = 16,
      "Worst" = 17,
      "Mid" = 18
    )
  ) +
  ylab("Annual hospitalisations per 100,000")+
  xlab("Age group")+
  theme_bw(base_size=14)+
  theme(panel.grid.minor = element_blank())

hosp_dat

ggsave(file.path("Figures", "hosp_data_best_worst.png"), height = 10, width = 14, dpi = 300)
  
###############################################################################################################
# Load the AIHW data 

df_hosp_aihw <- read_csv(file.path("data", "all_aihw_hosps_raw_2023.csv")) %>%
  filter(Age == "All") %>%
  rename(age_group = Age, hospitalisation = Separations,  hpht = Rate)


###############################################################################################################
# Save table to allow calculation of multiplier from FluCan to AIHW

df_hosp_flucan_filtered <- hosp_aggregated_flucan %>%
  filter(season == "2023") %>%
  select(-hpht) %>%
  group_by(season) %>%
  summarise(
    total_hospitalisations = sum(total_hosp),
    total_population = sum(total_pop),
    .groups = "drop"
  ) %>%
  mutate(
    hpht_flucan = total_hospitalisations * 100000 / total_population
  )

hosp_data_comparison_tbl <- df_hosp_aihw %>%
  rename(
    total_pop_aihw = Population,
    total_hosp_aihw = hospitalisation,
    hpht_aihw = hpht
    ) %>%
  mutate(
    season = "2023",
    hpht_flucan = df_hosp_flucan_filtered$hpht_flucan
    ) %>%
  select(season, age_group, hpht_flucan, hpht_aihw)


write.csv(hosp_data_comparison_tbl, file.path("data", "hosp_data_comparison_tbl.csv"), row.names = FALSE)

mult_flucan_to_aihw <- hosp_data_comparison_tbl$hpht_aihw / hosp_data_comparison_tbl$hpht_flucan

best_hosp_flucan_final <- best_hosp_flucan %>%
  mutate(
    hpht_flucan_best = hpht * mult_flucan_to_aihw,
    age_index = case_when(
        age == "0-1" ~ 1,
        age == "1-4" ~ 2,
        age == "5-11" ~ 3,
        age == "12-17" ~ 4,
        age == "18-64" ~ 5,
        age == "65-79" ~ 6,
        age == ">80" ~ 7
    )
  ) %>%
  select(hpht_flucan_best, age_index)
  
worst_hosp_flucan_final <- worst_hosp_flucan %>%
  mutate(
    hpht_flucan_worst = hpht * mult_flucan_to_aihw,
    age_index = case_when(
      age == "0-1" ~ 1,
      age == "1-4" ~ 2,
      age == "5-11" ~ 3,
      age == "12-17" ~ 4,
      age == "18-64" ~ 5,
      age == "65-79" ~ 6,
      age == ">80" ~ 7
    )
  ) %>%
  select(hpht_flucan_worst, age_index)
  
write.csv(best_hosp_flucan_final, file.path("data", "best_hosp_flucan_final.csv"), row.names = FALSE)  
write.csv(worst_hosp_flucan_final, file.path("data", "worst_hosp_flucan_final.csv"), row.names = FALSE)  





