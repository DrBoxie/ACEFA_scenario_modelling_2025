################################################################################.
#                               Colour palettes                                # -----
################################################################################.
#install.packages("devtools")
#devtools::install_github("idem-lab/idpalette")
library(idpalette) # https://github.com/acefa-hubs/idpalette
library(colorspace)

# idpal('acefa')[c(1,3,5)]
# idpalette("acefa", 9)[8] # "#E375BD"


scenario_palettes <- list()
scenario_palettes$team <- c("UoM" = "#262261",
                            "UNSW" = "#eb008b")
scenario_palettes$team <- c("no waning" = "#262261",
                            "with waning" = "#eb008b")
scenario_palettes$setting <- c("temperate" = "#262261",
                               "tropical" = "#eb008b")

scenario_palettes$scenario <- c("Status_quo" = "#662d91", 
                                "Pessimistic_LAIV_5-12yo" = "#84b5cc", 
                                "Mid_LAIV_5-12yo" = "#58b4ac",
                                "Optimistic_LAIV_5-12yo" = "#E375BD",
                                
                                "Pessimistic_LAIV_5-18yo" = "#262261",
                                "Mid_LAIV_5-18yo" = "#178793",
                                "Optimistic_LAIV_5-18yo" = "#eb008b",
                                "Mid_LAIV_2-5yo" = "#4aa3a0",  
                                "Mid_LAIV_2-12yo" = "#2e8f8c",  
                                "Mid_LAIV_2-18yo" = "#1f7c78",
                                
                                "Mid_Low_LAIV_5-12yo" = "#3b9a97",
                                "Mid_Low_LAIV_5-18yo" = "#7bc9c2",
                                "Mid_Low_LAIV_2-12yo" = "#6cc0b8",
                                "Mid_Low_LAIV_2-18yo" = "#16706d",
                                "Mid_Extra_Low_LAIV_5-12yo" = "#135f5c",
                                "Mid_Extra_Low_LAIV_5-18yo" = "#0f4e4c"
                                )

scenario_palettes$agegroup <- c("All" = "black",
                                "<1" = "#eb008b",
                                "1-4" = "#E375BD",
                                "5-11" = "#84b5cc",
                                "12-17" = "#58b4ac",
                                "18-64"= "#178793",
                                "65-79" = "#262261",
                                "80+"= "#662d91")


labeller_target <- c("infection_incidence" = "Daily infection incidence",
                     "disease_incidence" = "Daily case incidence",
                     "admission_incidence" = "Daily hospital admissions")


labeller_scenario <- c("Status_quo" = "Status quo", 
                       "Pessimistic_LAIV_5-12yo" = "LAIV ages 5-<12\n(Pessimistic effectiveness)", 
                       "Mid_LAIV_5-12yo" = "LAIV ages 5-<12\n(Central effectiveness)",
                       "Optimistic_LAIV_5-12yo" = "LAIV ages 5-<12\n(Optimistic effectiveness)",
                       "Pessimistic_LAIV_5-18yo" = "LAIV ages 5-<18\n(Pessimistic effectiveness)",
                       "Mid_LAIV_5-18yo" = "LAIV ages 5-<18\n(Central effectiveness)",
                       "Optimistic_LAIV_5-18yo" = "LAIV ages 5-<18\n(Optimistic effectiveness)",
                       "Mid_LAIV_2-5yo" = "LAIV ages 2-<5\n(Central effectiveness)",
                       "Mid_LAIV_2-12yo" = "LAIV ages 2-<12\n(Central effectiveness)",
                       "Mid_LAIV_2-18yo" = "LAIV ages 2-<18\n(Central effectiveness)",
                       
                       "Mid_Low_LAIV_5-12yo" = "LAIV ages 5-<12\n(Central effectiveness)\n(Low coverage)",
                       "Mid_Low_LAIV_5-18yo" = "LAIV ages 5-<18\n(Central effectiveness)\n(Low coverage)",
                       "Mid_Low_LAIV_2-12yo" = "LAIV ages 2-<12\n(Central effectiveness)\n(Low coverage)",
                       "Mid_Low_LAIV_2-18yo" = "LAIV ages 2-<18\n(Central effectiveness)\n(Low coverage)",
                       "Mid_Extra_Low_LAIV_5-12yo" = "LAIV ages 5-<12\n(Central effectiveness)\n(Extra low coverage)",
                       "Mid_Extra_Low_LAIV_5-18yo" = "LAIV ages 5-<18\n(Central effectiveness)\n(Extra low coverage)"
                       )

order_scenario <- c("Status_quo","Pessimistic_LAIV_5-12yo", "Mid_LAIV_5-12yo",
                    "Optimistic_LAIV_5-12yo","Pessimistic_LAIV_5-18yo" ,"Mid_LAIV_5-18yo","Optimistic_LAIV_5-18yo",
                    "Mid_LAIV_2-5yo", "Mid_LAIV_2-12yo", "Mid_LAIV_2-18yo", "Mid_Low_LAIV_5-12yo",  "Mid_Low_LAIV_5-18yo",
                    "Mid_Low_LAIV_2-12yo", "Mid_Low_LAIV_2-18yo", "Mid_Extra_Low_LAIV_5-12yo",  "Mid_Extra_Low_LAIV_5-18yo" )

# Add in blank panels
order_scenario2 <- c("", "Status_quo", " ", "Pessimistic_LAIV_5-12yo", "Mid_LAIV_5-12yo",
                     "Optimistic_LAIV_5-12yo","Pessimistic_LAIV_5-18yo" ,"Mid_LAIV_5-18yo","Optimistic_LAIV_5-18yo",
                     "Mid_LAIV_2-5yo", "Mid_LAIV_2-12yo", "Mid_LAIV_2-18yo", "Mid_Low_LAIV_5-12yo",  "Mid_Low_LAIV_5-18yo",
                     "Mid_Low_LAIV_2-12yo", "Mid_Low_LAIV_2-18yo", "Mid_Extra_Low_LAIV_5-12yo",  "Mid_Extra_Low_LAIV_5-18yo" )

darker_team_palette <- sapply(scenario_palettes$team, darken, 0.3)

