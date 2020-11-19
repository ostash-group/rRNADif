# Pased arguments
# $1 - The file which to use for plot
# $2 - The working directory

library(tidyr)
library(dplyr)
library(ggplot2)

# Get all the command line arguments
args = commandArgs(trailingOnly=TRUE)
# Set working directory
setwd(args[2])

# Read the data
df_plot <- read.csv(args[1])

# Plot the data
plot <- df_plot %>%
  pivot_longer(c("Mean", "Median"),names_to = "type", values_to = "values") %>%
  ggplot(aes(x=values, group=type, fill=type)) +
  geom_density(adjust=1.5, alpha=.4)+
  xlab("Mean and Median branch length values") +
  ylab("Density")
  
# Save the plot
ggsave("Plot_for_non_outliers.pdf", plot = plot,width = 7 , height = 5 ,dpi = 1200)