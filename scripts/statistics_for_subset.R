#Passed arguments:
# $1 - The database name for proper filenames of the results
# $2 - Generated subset of precomputed database (from 
#      dataframe_preprocess_for_subset.py)
# $3 - Directory from which the script should work (here is tmp_database folder)
# $4 - The folder where to write result files

library(ape)

# Get the arguments
args = commandArgs(trailingOnly=TRUE)

# Set the script execution directory:
setwd(args[3])

# Read the input dataframe (subset from master)
df_mean_median = read.csv(args[2])

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
write.csv(df_mean_median,paste0(args[4],"/",args[1], "_all.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0(args[4],"/",args[1], "_median_outliers.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0(args[4],"/",args[1], "_mean_outliers.csv"), row.names = FALSE)
write.csv(df_norm,paste0(args[4],"/",args[1], "_no_ouliers.csv"), row.names = FALSE)
