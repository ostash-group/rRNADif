# Passed arguments:
# $1 - dataframe for the input
# $2 - Database dataframe
# $3 - Directory were the files are present
library(ape)

# Get all the arguments
args = commandArgs(trailingOnly=TRUE)

# Setting working directory
setwd(args[3])

# Read dataframes and combine them
df_1 <- read.csv(args[1])
df_2 <- read.csv(args[2])
df_mean_median <- rbind(df_1, df_2)

# Get the mean/median outliers 
outliers_mean <- boxplot(df_mean_median$Mean)$out
out_mean <- unique(df_mean_median[which(df_mean_median$Mean %in% outliers_mean),]$Names)
outliers_median <- boxplot(df_mean_median$Median)$out
out_median <- unique(df_mean_median[which(df_mean_median$Median %in% outliers_median),]$Names)

# Make dataframes out of them 
df_mean_outliers <- df_mean_median[df_mean_median$Names %in% out_mean,]
df_median_outliers <- df_mean_median[df_mean_median$Names %in% out_median,]

# Filter the dataframe to remain inly non-outlier organisms
df_norm <- df_mean_median[!(df_mean_median$Names %in% out_median),]
df_norm <- df_norm[!(df_norm$Names %in% out_mean),]

# Write results into dataframes
write.csv(df_mean_median,paste0("Results_all.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0("Results_median_outliers.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0("Results_mean_outliers.csv"), row.names = FALSE)
write.csv(df_norm,paste0("Results_no_ouliers.csv"), row.names = FALSE)