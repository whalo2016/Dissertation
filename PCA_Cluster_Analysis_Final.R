library(readr)
library(dplyr)
library(factoextra)
library(mclust)
library(purrr)
library(ggplot2)
library(cluster)

#____________________________Part 1 - PCA Cluster Analysis_________________________________
#//////////////////////////////////////////////////////////////////////////////////////////

#Open Attributes File
CatchmentAttributesNSE <- read.csv("CatchmentAttributesNSE.csv", skip = 1, header = TRUE)

#------------------------------------------------------------------------------------------
#Section 1 - Data Cleaning and Scaling

# 1. Remove First 7 Columns
attribute_data <- CatchmentAttributesNSE %>%
  select(-c(1,2,3,4,5,6,7))

# 2. Keep only numeric columns
attribute_data <- attribute_data[sapply(attribute_data, is.numeric)]

# 3. Remove columns with all NA or constant values
attribute_data <- attribute_data[, colSums(is.na(attribute_data)) < nrow(attribute_data)] 
attribute_data <- attribute_data[, apply(attribute_data, 2, sd, na.rm = TRUE) != 0]     

# 4. Replace remaining NA/Inf with column means
for(i in seq_along(attribute_data)) {
  col <- attribute_data[, i]
  col[!is.finite(col)] <- NA
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  attribute_data[, i] <- col
}

# 5. Scale the cleaned data
scaled_data <- scale(attribute_data)

#------------------------------------------------------------------------------------------
#Section 2 - PCA

# 6. PCA Analysis
pca_results <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
loadings <- pca_results$rotation
write.csv(loadings, "PCLoadings.csv")

var_explained <- (pca_results$sdev^2) / sum(pca_results$sdev^2)
cumulative_var <- cumsum(var_explained)

which(cumulative_var >= 0.90)[1]

pca_table <- data.frame( 
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained,
  Cumulative = cumulative_var ) 
pca_table

# 7. Scree Plot
fviz_eig(pca_results, addlabels = TRUE, barfill = "lightblue", barcolor = "lightblue") +
  labs(
    title = "Scree Plot of Principal Components",
    x = "Principal Components",
    y = "Percentage of Explained Variance"
  ) +
  theme_classic(base_size = 14)

# 8. Keep first 10 PCs for clustering(based on scree plot and cumulative variance)
pc_data <- pca_results$x[, 1:10]

# 9. Identify variables with greatest positive and negative loadings for PC1-10
for (pc in colnames(loadings)[1:10]) {
  cat("\n", pc, "\n")
  
  TopLoadings <- data.frame(var = rownames(loadings),
                            loading = loadings[, pc],
                            row.names = NULL)
  
  # Top 5 positive
  cat("Top positive loadings:\n")
  print(head(TopLoadings[order(-TopLoadings$loading), ], 5))
  
  # Top 5 negative
  cat("Top negative loadings:\n")
  print(head(TopLoadings[order(TopLoadings$loading), ], 5))
}

# 10. PCA Variable Plot
fviz_pca_var(pca_results, col.var = "contrib")

# 11. PCA Individual Plot
fviz_pca_ind(pca_results,
             geom.ind = "point",
             col.ind = "blue",
             repel = TRUE,
             title = "Gauges plotted on PC1 vs PC2")

#------------------------------------------------------------------------------------------
# Section 3 - Determining Optimal K

# 12. Elbow Method to identify optimal k
set.seed(123)
fviz_nbclust(pc_data, kmeans, method = "wss") +
  labs(title = "Elbow Method for Optimal k (PCA Clustering)",
       x = "Number of clusters (k)",
       y = "Total Within-Cluster Sum of Squares") +
  theme_classic(base_size = 14)

# 13. Silhouette method
fviz_nbclust(pc_data, kmeans, method = "silhouette") +
  labs(title = "Silhouette Analysis for Optimal k (PCA Clustering)",
       x = "Number of Clusters (k)",
       y = "Average Silhouette Width") +
  theme_classic(base_size = 14)

# 14. Stability check with Adjusted Rand Index (ARI)
set.seed(123)
k_range <- 2:10
ari_scores <- sapply(k_range, function(k){
  km1 <- kmeans(scaled_data, centers = k, nstart = 50)
  km2 <- kmeans(scaled_data, centers = k, nstart = 50)
  adjustedRandIndex(km1$cluster, km2$cluster)
})

# 15. Plot Stability Check with ARI
par(font.main = 1,   
    font.lab = 1,    
    cex.main = 1.3,  
    cex.lab = 1.2,  
    cex.axis = 1.1)  

plot(k_range, ari_scores, type = "b",
     pch = 19,
     col = "lightblue",
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index (stability)",
     main = "Cluster Stability Across Repeated Runs (PCA Clustering)")


#------------------------------------------------------------------------------------------
# Section 4 - K-means Clustering

# 16. Run k-means with optimal k on reduced PC data
set.seed(123)
km_pca <- kmeans(pc_data, centers = 3, nstart = 25)

# 17. Attach cluster labels
attribute_data$cluster <- km_pca$cluster

clusters <- km_pca$cluster
table(clusters)

# 18. PCA plot coloured by cluster
fviz_pca_ind(pca_results,
             geom.ind = "point",
             col.ind = as.factor(clusters),   # color by cluster
             palette = "jco",                 # nice color palette
             addEllipses = TRUE,              # draw cluster ellipses
             legend.title = "Cluster") +
  labs(
    title = "Gauges in PCA space colored by clusters",
    x = "PC1 (32.5%)",
    y = "PC2 (17.3%)") +
  theme_classic()

# 19. View cluster centers
km_pca$centers

# 20. Count how many catchments per cluster
table(clusters)

# 21. Merge clusters back to full dataset
CatchmentAttributesNSE$cluster <- km_pca$cluster

# 22. Export selected performance metrics and cluster labels
export_data <- CatchmentAttributesNSE[, c("gauge_id", "NSE", "RMSE", "R2", "n_points", "cluster")]
write.csv(export_data, "Cluster_Performance.csv")

#------------------------------------------------------------------------------------------
#Section 5 - Statistical Correlation Tests: NSE vs PC1-10

# 23. Add NSE to PC dataframe
pc_df <- as.data.frame(pc_data)
pc_df$NSE <- CatchmentAttributesNSE$NSE

# 24. Spearman Correlations for PC1-10
pc_spearman <- map_df(names(pc_df)[1:10], function(pc) {
  
  x <- pull(pc_df, pc)
  y <- pull(pc_df, NSE)
  
  broom::tidy(cor.test(x, y, method = "spearman")) %>%
    mutate(PC = pc)
})
pc_spearman


# 25. Linear Models NSE - PC Score
pc_linear <- map_df(names(pc_df)[1:10], function(pc) { 
  broom::tidy(lm(NSE ~ pc_df[[pc]], data = pc_df)) %>% 
    mutate(PC = pc)
})
pc_linear

#_____________________________Part 2 - Cluster Performance_________________________________
#//////////////////////////////////////////////////////////////////////////////////////////

# 26. Open Performance File 
performance_data <- read.csv("Cluster_Performance.csv")

#------------------------------------------------------------------------------------------
#Section 1 - Performance Summary Cluster-by-Cluster

# 27. Summarise model performance by cluster
performance_summary <- performance_data %>%
  group_by(cluster) %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n()
  )
print(performance_summary)

#------------------------------------------------------------------------------------------
#Section 2 - Visualise NSE Distributions

# 28. Visualise NSE distributions per cluster
boxplot(NSE ~ cluster, data = performance_data,
        main = "NSE distribution by cluster",
        xlab = "Cluster",
        ylab = "NSE",
        col = "lightblue")


ggplot(performance_data, aes(x = factor(cluster), y = NSE, fill = factor(cluster))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "NSE distribution per cluster",
       x = "Cluster", y = "NSE") +
  theme_minimal()

#------------------------------------------------------------------------------------------
#Section 3 - Outlier Extraction - Data Filtering

#Can clearly see massive outlier in cluster 2 and 3
#Identify what gauge has this outlier
performance_data %>% filter(cluster == 2, NSE < -1)
performance_data %>% filter(cluster == 3, NSE < -1)


#Redo performance summary after removing outlier
filtered_data <- performance_data  %>% filter(gauge_id != 3050 & gauge_id != 3024)
performance_summary <- filtered_data %>%
  group_by(cluster) %>%
  summarise(
    mean_NSE = mean(NSE, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_NSE = sd(NSE, na.rm = TRUE),
    iqr_NSE = IQR(NSE, na.rm = TRUE),
    n_gauges = n()
  )
print(performance_summary)

#------------------------------------------------------------------------------------------
#Section 4 - Visualise NSE Distributions without outliers

#Visualise NSE distributions per cluster
boxplot(NSE ~ cluster, data = filtered_data,
        main = "Distribution of NSE Across Clusters",
        xlab = "Cluster",
        ylab = "Nash-Sutcliffe Efficiency (NSE)",
        col = "lightblue",
        border = "#2c3e50",
        outline = FALSE,
        cex.lab = 1.2, 
        cex.main = 1.3, 
        cex.axis = 1.1)
stripchart(NSE ~ cluster, data = filtered_data,
           method = "jitter",
           pch = 16,
           col = rgb(0,0,0,0.4),
           vertical = TRUE,
           add = TRUE)

#Create violin plot with mean and median overlays
ggplot(filtered_data, aes(x = factor(cluster), y = NSE, fill = factor(cluster))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Mean (black dot) and median (red triangle)
  geom_point(data = performance_summary, aes(x = factor(cluster), y = mean_NSE),
             color = "black", size = 3, shape = 16) +
  geom_point(data = performance_summary, aes(x = factor(cluster), y = median_NSE),
             color = "red", size = 3, shape = 17) +
  # Labels for numeric values
  geom_text(data = performance_summary,
            aes(x = factor(cluster), y = mean_NSE, label = round(mean_NSE, 2)),
            vjust = 1.5, color = "black", size = 3) +
  geom_text(data = performance_summary,
            aes(x = factor(cluster), y = median_NSE, label = round(median_NSE, 2)),
            vjust = -1.5, color = "red", size = 3) +
  labs(title = "NSE distribution per cluster (outliers removed)",
       x = "Cluster", y = "NSE") +
  theme_minimal()



ggplot(filtered_data, aes(x = factor(cluster), y = NSE, fill = factor(cluster))) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Mean (black dot) and median (red triangle) with legend mapping
  geom_point(data = performance_summary, aes(x = factor(cluster), y = mean_NSE, colour = "Mean", shape = "Mean"),
    size = 3) +
  geom_point(data = performance_summary, aes(x = factor(cluster), y = median_NSE, colour = "Median", shape = "Median"),
    size = 3) +
  
  # Labels for numeric values
  geom_text(data = performance_summary,
    aes(x = factor(cluster), y = mean_NSE, label = round(mean_NSE, 2)),
    vjust = 1.5, colour = "black", size = 3) +
  geom_text(data = performance_summary,
    aes(x = factor(cluster), y = median_NSE, label = round(median_NSE, 2)),
    vjust = -1.5, colour = "red", size = 3) +
  
  labs(title = "Distribution of NSE Across Clusters (outliers removed)",
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


#------------------------------------------------------------------------------------------
#Section 5 - Statistical Tests (DO clusters statistically differ in NSE?)

#ANOVA
anova_result <- aov(NSE ~ factor(cluster), data = filtered_data)
summary(anova_result)

# Kruskal–Wallis (non-parametric, safer for NSE) 
kruskal_result <- kruskal.test(NSE ~ factor(cluster), data = filtered_data) 
kruskal_result







