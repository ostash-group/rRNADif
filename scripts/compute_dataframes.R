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
setwd(args[2])

files <- list.files( pattern="*.nwk$", full.names=TRUE, recursive=FALSE)
a <- list()
mean_median <- list()


for(i in 1:length(files)){
  tree <- read.tree(files[i])
  matrix <- try(cophenetic.phylo(tree))
  if ("try-error" %in% class(matrix)) {
    i <-  i+1
    tree <- read.tree(files[i])
    matrix <- cophenetic.phylo(tree)
  }
  col_name <- rownames(matrix)[which(nchar(rownames(matrix)) == min(nchar(rownames(matrix))))]
  matrix[lower.tri(matrix)] <- NA
  diag(matrix) <- NA
  m_vector <- c(matrix)
  m_vector<-m_vector[!is.na(m_vector)]
  a[[col_name]] <- m_vector
  mean_median[[col_name]] <- c(mean(m_vector), median(m_vector))
}


Names <- c()
Values <- c()
for (i in 1:length(a)){
  Names <- c(Names, rep(names(a[i]), length(a[[i]])))
  Values <- c(Values, a[[i]])
}

Names_1 <- c()
mean_1 <- c()
median_1 <- c()
for (i in 1:length(a)){
  Names_1 <- c(Names_1, names(mean_median[i]))
  mean_1 <- c(mean_1, mean_median[[i]][1])
  median_1 <- c(median_1, mean_median[[i]][2])
}


df <- data.frame(Names, Values)
df_mean_median <- data.frame(Names_1, mean_1, median_1)

colnames(df_mean_median) <- c("Names", "Mean", "Median")

outliers_mean <- boxplot(df_mean_median$Mean)$out
out_mean <- unique(df_mean_median[which(df_mean_median$Mean %in% outliers_mean),]$Names)

outliers_median <- boxplot(df_mean_median$Median)$out
out_median <- unique(df_mean_median[which(df_mean_median$Median %in% outliers_median),]$Names)

df_mean_outliers <- df_mean_median[df_mean_median$Names %in% out_mean,]

df_median_outliers <- df_mean_median[df_mean_median$Names %in% out_median,]

df_norm <- df_mean_median[!(df_mean_median$Names %in% out_median),]
df_norm <- df_norm[!(df_norm$Names %in% out_mean),]

write.csv(df_mean_median,paste0(args[1], "_all.csv"), row.names = FALSE)
write.csv(df_median_outliers,paste0(args[1], "_median_outliers.csv"), row.names = FALSE)
write.csv(df_mean_outliers,paste0(args[1], "_mean_outliers.csv"), row.names = FALSE)
write.csv(df_norm,paste0(args[1], "_no_ouliers.csv"), row.names = FALSE)
