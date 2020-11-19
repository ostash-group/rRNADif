# Passed arguments:
# $1 - Filename for proper result files naming
# $2 - Working directory

library(ape)

# Gather all the arguments from command line
args = commandArgs(trailingOnly=TRUE)

# Setting the workind directory:
setwd(args[2])

# List all *.nwk files
files <- list.files( pattern="*.nwk$", full.names=TRUE, recursive=FALSE)
# Create emtpty lists
a <- list()
mean_median <- list()

# Loop through files and extract branch lenght. Compute mean and median for 
# values and put into list. Use or organism name as name in a list
for(i in 1:length(files)){
  # REad the file
  tree <- read.tree(files[i])
  # Extract branch lengths
  matrix <- cophenetic.phylo(tree)
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


# Make vectors out of a list. One for names and two for values 
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

df_mean_median$Names <- args[1]

write.csv(df_mean_median,paste0(args[1], "_all.csv"), row.names = FALSE)