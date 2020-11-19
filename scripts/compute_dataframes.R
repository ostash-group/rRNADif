# Passed agruments:
# $1 - Name of a database. This will be used for file naming in the end
# $2 - The directory script should work in 

library(ape)

# Get all the arguments from commandline
args = commandArgs(trailingOnly=TRUE)

# Set the working directory 
setwd(args[2])

# Get list of files with .nwk extension in this directory
files <- list.files( pattern="*.nwk$", full.names=TRUE, recursive=FALSE)
# Create emtpty lists
a <- list()
mean_median <- list()

# Loop through files and extract branch lenght. Compute mean and median for 
# values and put into list. Use or organism name as name in a list
for(i in 1:length(files)){
  # Read the file
  tree <- read.tree(files[i])
  # Extract branch lengths
  matrix <- try(cophenetic.phylo(tree))
  # If smth wrong with the file jump to the next one (TO-DO: more reliable 
  # way of hangling this error)
  if ("try-error" %in% class(matrix)) {
    i <-  i+1
    tree <- read.tree(files[i])
    matrix <- cophenetic.phylo(tree)
  }
  # col_name is actually organism name (the shortest from all sequences' names,
  # as this would not have numbers after it(the numbers are used to omit
  # duplicate sequence names))
  col_name <- rownames(matrix)[which(nchar(rownames(matrix)) == min(nchar(rownames(matrix))))]
  # Set NA for diagonal and one of the triangles of a matrix. (Drop those values)
  matrix[lower.tri(matrix)] <- NA
  diag(matrix) <- NA
  # Compute mean and median from a vector
  m_vector <- c(matrix)
  m_vector<-m_vector[!is.na(m_vector)]
  a[[col_name]] <- m_vector
  mean_median[[col_name]] <- c(mean(m_vector), median(m_vector))
}

Names_1 <- c()
mean_1 <- c()
median_1 <- c()
for (i in 1:length(a)){
  Names_1 <- c(Names_1, names(mean_median[i]))
  mean_1 <- c(mean_1, mean_median[[i]][1])
  median_1 <- c(median_1, mean_median[[i]][2])
}

df_mean_median <- data.frame(Names_1, mean_1, median_1)

colnames(df_mean_median) <- c("Names", "Mean", "Median")

# Get the outlier mean/median statistics for a dataframe and put
# them into dataframe
outliers_mean <- boxplot(df_mean_median$Mean)$out
out_mean <- unique(df_mean_median[which(df_mean_median$Mean %in% outliers_mean),]$Names)
outliers_median <- boxplot(df_mean_median$Median)$out
out_median <- unique(df_mean_median[which(df_mean_median$Median %in% outliers_median),]$Names)
df_mean_outliers <- df_mean_median[df_mean_median$Names %in% out_mean,]
df_median_outliers <- df_mean_median[df_mean_median$Names %in% out_median,]

# Store non-outlier values in a different dataframe
df_norm <- df_mean_median[!(df_mean_median$Names %in% out_median),]
df_norm <- df_norm[!(df_norm$Names %in% out_mean),]

# Write the results
write.csv(df_mean_median,paste0(args[1], "_all.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0(args[1], "_median_outliers.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0(args[1], "_mean_outliers.csv"), row.names = FALSE)
write.csv(df_norm,paste0(args[1], "_no_ouliers.csv"), row.names = FALSE)
