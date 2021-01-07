#Passed arguments:
# $1 - Dataframe to compute
# $2 - Directory from which the script should work (here is tmp_database folder)
# $3 - Database name

library(ape)

# Get the arguments
args = commandArgs(trailingOnly=TRUE)

# Set the script execution directory:
setwd(args[2])

# Read the input dataframe (subset from master)
df_mean_median = read.csv(args[1])

# Generate outlier statistics fore means and median
outliers_mean <- boxplot(df_mean_median$Mean)$out
out_mean <- unique(df_mean_median[which(df_mean_median$Mean %in% outliers_mean),]$Names)
outliers_median <- boxplot(df_mean_median$Median)$out
out_median <- unique(df_mean_median[which(df_mean_median$Median %in% outliers_median),]$Names)

# Make smaller dataframes with mean and median outliers
df_mean_outliers <- df_mean_median[df_mean_median$Names %in% out_mean,]
df_median_outliers <- df_mean_median[df_mean_median$Names %in% out_median,]

# Make a dataframe with no outliers
df_norm <- df_mean_median[!(df_mean_median$Names %in% out_median),]
df_norm <- df_norm[!(df_norm$Names %in% out_mean),]

# Write result files
write.csv(df_mean_median,paste0(args[3], "_all_new.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0(args[3], "_median_outliers_new.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0(args[3], "_mean_outliers_new.csv"), row.names = FALSE)
write.csv(df_norm,paste0(args[3], "_no_ouliers_new.csv"), row.names = FALSE)
