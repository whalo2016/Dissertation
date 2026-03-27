library(dplyr)
library(readr)
library(purrr)
library(lubridate)
library(tidyr)

#Outlining file path for each catchment file for forcings

path <- "CAMELS_IND_All_Catchments/catchment_mean_forcings/"
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
length(files)

#Identifying the columns to keep from each file for model

cols_to_keep <- c("year",
                  "month",
                  "day",
                  "prcp(mm/day)",
                  "tavg(C)",
                  "aet_gleam(mm/day)"
                  )

#Iterative function, goes through each of the identified files
#Returns all the data needed from each file

cleaning <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  
  df <- df %>%
    select(any_of(cols_to_keep)) %>%
    mutate(
      date = make_date(year, month, day),
      catchment_id = tools::file_path_sans_ext(basename(path)),
      catchment_id = as.integer(catchment_id)) %>%
    select(-year, -month, -day) %>%
    select(date, everything())
           
  return(df)
}

#Adds all the data to a data frame and then writes to CSV
godavari_all <- map_df(files, cleaning)

write.csv(godavari_all, "Godavari_Forcings.csv")


#Streamflow CSV
df <- read_csv("streamflow_observed.csv", show_col_types = FALSE)

df <- df %>%
  mutate(date = make_date(year, month, day)) %>%
  select(-year, -month, -day) %>%
  select(date, everything())

df_formatted <- df %>%
  pivot_longer(cols = -date, names_to = "catchment_id", values_to = "streamflow")

df_formatted <- df_formatted %>%
  mutate(catchment_id = as.integer(catchment_id)) %>%
  filter(catchment_id >= 3001 & catchment_id <= 3110)

write_csv(df_formatted, "Godavari_Streamflow.csv")


#Joining both datasets together
forcings <- read_csv("Godavari_Forcings.csv", show_col_types = FALSE)
streamflow <- read_csv("Godavari_Streamflow.csv", show_col_types = FALSE)

#Cleaning Datasets before merging
forcings <- forcings %>%
  select(-1) %>%
  rename(
    prcp = `prcp(mm/day)`,
    tavg = `tavg(C)`,
    aet = `aet_gleam(mm/day)`
  ) %>%
  select(catchment_id, date, everything())

Godavari_combined <- forcings %>%
  left_join(streamflow, by = c("catchment_id", "date"))

write_csv(Godavari_combined, "Godavari_FinalData.csv")


head(Godavari_combined)
Godavari_combined %>%
  filter(catchment_id == "3110")





