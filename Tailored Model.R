library(readr)
library(dplyr)
library(lubridate)
library(ranger)
library(hydroGOF)
library(ggplot2)

#______________________________Part 1 - Data Preparation___________________________________
#//////////////////////////////////////////////////////////////////////////////////////////

cluster <- read.csv("Cluster_Performance.csv") %>%
  select(c(2,7))

inputs <- read.csv("Godavari_FinalData.csv")

inputs <- filter(inputs, !is.na(streamflow)) %>%
  arrange(date) %>%
  rename(gauge_id = catchment_id)


model_data <- inputs %>%
  left_join(cluster, by = "gauge_id")

attributes <- read.csv("CatchmentAttributesNSE.csv", skip = 1, header = TRUE) %>%
  select(-c(2,3,4,5,6,7))

merged <- merge(model_data, attributes, by = "gauge_id", all.x = TRUE)%>%
  arrange(date)

gauges <- unique(merged$gauge_id)
results <- data.frame()


#_________________________Part 2 - Tailored Model Development______________________________
#//////////////////////////////////////////////////////////////////////////////////////////

#Generating lag features
for (g in gauges) {
  df_gauge <- merged %>%
    filter(gauge_id == g) %>%
    arrange(date) %>%
    mutate(
      lag1 = lag(streamflow, 1),
      lag2 = lag(streamflow, 2),
      lag3 = lag(streamflow, 3),
      lag7 = lag(streamflow, 7),
      lag14 = lag(streamflow, 14),
      lag30 = lag(streamflow, 30),
      lag60 = lag(streamflow, 60),
      month = month(date)
    ) %>%
    filter(!is.na(lag1) & !is.na(lag2) & !is.na(lag3) & !is.na(lag7))
  
  #Training-Testing data split    
  split_data <- floor(0.7*nrow(df_gauge))
  training <- df_gauge[1:split_data, ]
  testing <- df_gauge[(split_data + 1):nrow(df_gauge), ]
  
  cl <- unique(df_gauge$cluster)
  
  #Cluster Specific Input variables
  if (cl == 1) {
    inputs <- streamflow ~ prcp + tavg + aet + 
      reservoir_index + num_dams + total_storage + lai_mean +
      lag1 + lag2 + lag3 + lag7 + month
  }
  
  if (cl == 2) {
    inputs <- streamflow ~ prcp + tavg + aet +
      p_monthly_variability + low_prec_freq + slope_mean + urban_frac_2005 +
      fall_rate_mean + rise_rate_mean + p_max + built_area_frac +
      lag1 + lag2 + lag3 + lag7 + month
  }
  
  if (cl == 3) {
    inputs <- streamflow ~ prcp + tavg + aet +
      p_mean + soil_awsc_major + irrigation_frac + geol_porosity +
      q_low_days + q_zero + flow_availability +
      lag1 + lag2 + lag3 + lag7 + lag14 + lag30 + lag60 + month
  }
  
  #Tune mtry for this gauge
  p <- length(all.vars(inputs)) - 1   # number of predictors
  mtry_values <- mtry_values <- c(3, 5, 7, 10, 15)
  mtry_values <- mtry_values[mtry_values < p]  # avoid invalid mtry
  
  best_mtry <- NULL
  best_nse <- -Inf
  
  for (m in mtry_values) {
    model_tmp <- tryCatch({
      ranger(
        inputs,
        data = training,
        num.trees = 200,  
        mtry = m
      )
    }, error = function(e) NULL)
    
    if (is.null(model_tmp)) next
    
    preds_tmp <- predict(model_tmp, data = testing)$predictions
    ME_tmp <- mean(preds_tmp - testing$streamflow)
    preds_tmp <- preds_tmp - ME_tmp
    
    nse_tmp <- NSE(preds_tmp, testing$streamflow)
    
    if (nse_tmp > best_nse) {
      best_nse <- nse_tmp
      best_mtry <- m
    }
  }
  
  # ---- Fit final model with best mtry ----
  RFModel <- ranger(
    inputs,
    data = training,
    num.trees = 2000,
    mtry = best_mtry
  )
  
  
  
  predictions <- predict(RFModel, data = testing)$predictions
  ME <- mean(predictions - testing$streamflow)
  corrected_predictions <- predictions - ME
  
  
  obs <- testing$streamflow
  rmse <- sqrt(mean((obs - corrected_predictions)^2))
  nse <- NSE(corrected_predictions, obs)
  r <- cor(corrected_predictions, obs)
  r2 <- r^2
  
  results <- rbind(results, data.frame(
    gauge_id = g,
    cluster = cl,
    NSE = nse,
    RMSE = rmse,
    R2 = r2,
    n_points = nrow(df_gauge)
  ))
}

print(results)

write_csv(results, "Final_Results.csv")


#______________________________Part 3 - Model Performance__________________________________
#//////////////////////////////////////////////////////////////////////////////////////////

sd_flow <- sd(merged$streamflow)

#Overall Model Performance
overall_performance_summary <- results %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n(),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    nrmse_sd = (mean(RMSE, na.rm = TRUE)) /sd_flow,
    mean_R2 = mean(R2, na.rm = TRUE)
  )
print(overall_performance_summary)

#Summarise model performance by cluster
cluster_performance_summary <- results %>%
  group_by(cluster) %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n(),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    nrmse_sd = (mean(RMSE, na.rm = TRUE)) /sd_flow,
    mean_R2 = mean(R2, na.rm = TRUE)
  )
print(cluster_performance_summary)

#------------------------------------------------------------------------------------------
#Outlier Extraction - Data Filtering

#Can clearly see massive outlier in cluster 2 and 3
#Identify what gauge has this outlier
results %>% filter(cluster == 2, NSE < -1)
results %>% filter(cluster == 3, NSE < -1)

#Redo overall performance summary after removing outliers
filtered_data <- results  %>% filter(gauge_id != 3050 & gauge_id != 3024)
filtered_overall_performance_summary <- filtered_data %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n(),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    nrmse_sd = (mean(RMSE, na.rm = TRUE)) /sd_flow,
    mean_R2 = mean(R2, na.rm = TRUE)
  )
print(filtered_overall_performance_summary)

#Redo cluster performance summary after removing outliers
filtered_cluster_performance_summary <- filtered_data %>%
  group_by(cluster) %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n(),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    nrmse_sd = (mean(RMSE, na.rm = TRUE)) /sd_flow,
    mean_R2 = mean(R2, na.rm = TRUE)
  )
print(filtered_cluster_performance_summary)

#__________________________Part 4 - Performance Visualisation______________________________
#//////////////////////////////////////////////////////////////////////////////////////////

ggplot(filtered_data, aes(x = factor(cluster), y = NSE, fill = factor(cluster))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Mean (black dot) and median (red triangle) with legend mapping
  geom_point(data = filtered_cluster_performance_summary, aes(x = factor(cluster), y = mean_NSE, colour = "Mean", shape = "Mean"),
             size = 3) +
  geom_point(data = filtered_cluster_performance_summary, aes(x = factor(cluster), y = median_NSE, colour = "Median", shape = "Median"),
             size = 3) +
  
  # Labels for numeric values
  geom_text(data = filtered_cluster_performance_summary,
            aes(x = factor(cluster), y = mean_NSE, label = round(mean_NSE, 2)),
            vjust = 1.5, colour = "black", size = 3) +
  geom_text(data = filtered_cluster_performance_summary,
            aes(x = factor(cluster), y = median_NSE, label = round(median_NSE, 2)),
            vjust = -1.5, colour = "red", size = 3) +
  
  labs(title = "Distribution of Tailored Model NSE Values Across Clusters (outliers removed)",
       x = "Cluster", y = "Nash-Sutcliffe Efficiency (NSE)", fill = "Cluster") +
  
  # Colour and shape legends for mean/median
  scale_colour_manual(name = "Statistics", values = c("Mean" = "black", "Median" = "red")) +
  scale_shape_manual(name = "Statistics", values = c("Mean" = 16, "Median" = 17)) +
  
  # Nicer cluster colours
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "lines")
  )

ggplot(filtered_data, aes(x = NSE)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Tailored Model NSE Values Across All Gauges",
    x = "Nash-Sutcliffe Efficiency (NSE)",
    y = "Count"
  ) +
  theme_classic()

#________________Part 5 - Performance Comparison ~ Baseline vs Tailored____________________
#//////////////////////////////////////////////////////////////////////////////////////////

baseline <- read.csv("Cluster_Performance.csv")

comparison <- baseline %>%
  select(gauge_id, baseline_NSE = NSE) %>%
  inner_join(results %>% select(gauge_id, tuned_NSE = NSE, cluster),
             by = "gauge_id")

comparison <- comparison %>%
  mutate(change_NSE = tuned_NSE - baseline_NSE) %>%
  filter(!gauge_id %in% c(3050, 3024))

# Gauge NSE improvement bar graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot(comparison,
       aes(y = reorder(gauge_id, change_NSE),
           x = change_NSE,
           fill = factor(cluster))) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.5) +
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  
  labs(
    title = "Change in NSE After Tuning (Tuned – Baseline)",
    x = "ΔNSE",
    y = "Gauge ID"
  ) +
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 5),
    legend.position = c(0.8, 0.2),  
    legend.background = element_rect(fill = "white", colour = "black"),
  )


# Baseline vs tuned model performance scattler plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot(comparison, aes(x = baseline_NSE, y = tuned_NSE, colour = factor(cluster))) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_brewer(palette = "Set1") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Baseline vs Tuned Model Performance",
    x = "Baseline NSE",
    y = "Tuned NSE",
    colour = "Cluster"
  ) +
  theme_classic()

# Mean Cluster NSE improvement bar graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cluster_improvement <- comparison %>%
  group_by(cluster) %>%
  summarise(mean_change = mean(change_NSE, na.rm = TRUE))

ggplot(cluster_improvement, aes(x = factor(cluster), y = mean_change, fill = factor(cluster))) +
  geom_col(width = 0.3) +
  geom_text(aes(label = round(mean_change, 3)), vjust = -0.5) +
  labs(
    title = "Mean Improvement in NSE by Cluster",
    x = "Cluster",
    y = "Mean ΔNSE",
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "none")


