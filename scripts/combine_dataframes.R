library(ape)


args = commandArgs(trailingOnly=TRUE)


getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

# Setting the script path would then be:
setwd(args[3])

df_1 <- read.csv(args[1])
df_2 <- read.csv(args[2])
df_mean_median <- rbind(df_1, df_2)

outliers_mean <- boxplot(df_mean_median$Mean)$out
out_mean <- unique(df_mean_median[which(df_mean_median$Mean %in% outliers_mean),]$Names)

outliers_median <- boxplot(df_mean_median$Median)$out
out_median <- unique(df_mean_median[which(df_mean_median$Median %in% outliers_median),]$Names)

df_mean_outliers <- df_mean_median[df_mean_median$Names %in% out_mean,]

df_median_outliers <- df_mean_median[df_mean_median$Names %in% out_median,]

df_norm <- df_mean_median[!(df_mean_median$Names %in% out_median),]
df_norm <- df_norm[!(df_norm$Names %in% out_mean),]

write.csv(df_mean_median,paste0("Results_all.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0("Results_median_outliers.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0("Results_mean_outliers.csv"), row.names = FALSE)
write.csv(df_norm,paste0("Results_no_ouliers.csv"), row.names = FALSE)
