library(dplyr)
library(lubridate)
library(hydroGOF)
library(ranger)
library(readr)

df <- read_csv("Godavari_FinalData.csv", show_col_types = FALSE)

df <- filter(df, !is.na(streamflow)) %>%
  arrange(date)

catchments <- unique(df$catchment_id)
results <- data.frame()

#_________________________Part 1 - Baseline Model Development______________________________
#//////////////////////////////////////////////////////////////////////////////////////////

#Generating lag features
for (cid in catchments) {
  df_catchment <- df %>% 
    filter(catchment_id == cid) %>%
    mutate(
      lag1 = lag(streamflow, 1),
      lag2 = lag(streamflow, 2),
      lag3 = lag(streamflow, 3),
      lag7 = lag(streamflow, 7),
      month = month(date)
    ) %>%
    filter(!is.na(lag1) & !is.na(lag2) & !is.na(lag3) & !is.na(lag7))

  #Training-Testing data split    
  split_data <- floor(0.7*nrow(df_catchment))
  training <- df_catchment[1:split_data, ]
  testing <- df_catchment[(split_data + 1):nrow(df_catchment), ]
  
  #Input of predictor variables into RF model
  #Hyperparameter selection
  set.seed(123)
  RFModel <- ranger(
    streamflow ~ prcp + tavg + aet + lag1 + lag2 + lag3 + lag7 + month,
    data = training,
    num.trees = 2000,
    mtry = 4
  )
  
  #Performance Metrics
  predictions <- predict(RFModel, data = testing)$predictions
  ME <- mean(predictions - testing$streamflow)
  corrected_predictions <- predictions - ME
  
  
  obs <- testing$streamflow
  rmse <- sqrt(mean((obs - corrected_predictions)^2))
  nse <- NSE(corrected_predictions, obs)
  r <- cor(corrected_predictions, obs)
  r2 <- r^2
  
  results <- rbind(results, data.frame(
    catchment_id = cid,
    NSE = nse,
    RMSE = rmse,
    R2 = r2,
    n_points = nrow(df_catchment)
  ))
}

results %>% arrange(desc(NSE))
print(results)

write_csv(results, "CatchmentNSE_Results.csv")


results <- read.csv("CatchmentNSE_Results.csv")


#Creation of boxplot
min_flow <- min(df$streamflow)
max_flow <- max(df$streamflow)
mean_flow <- mean(df$streamflow)
sd_flow <- sd(df$streamflow)

nrmse <- 536.704 / (max_flow - min_flow)
nrmse_mean <- 536.704 / mean_flow
nrmse_sd <- 536.704 / sd_flow

nrmse_sd_outlier <- 553.750 / sd_flow

nrmse_mean
nrmse_sd
nrmse_sd_outlier

boxplot(results$NSE,
        ylab = "Nash-Sutcliffe Efficiency (NSE)",
        main = "Distribution of NSE Across Catchments",
        col = "lightblue",
        border = "#2c3e50",
        outline = FALSE,
        cex.lab = 1.2, 
        cex.main = 1.3, 
        cex.axis = 1.1)

stripchart(results$NSE,
           method = "jitter",
           pch = 16,
           col = rgb(0,0,0,0.4),
           vertical = TRUE,
           add = TRUE)
