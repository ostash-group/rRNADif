library(tidyr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[2])

df_plot <- read.csv(args[1])

plot <- df_plot %>%
  pivot_longer(c("Mean", "Median"),names_to = "type", values_to = "values") %>%
  ggplot(aes(x=values, group=type, fill=type)) +
  geom_density(adjust=1.5, alpha=.4)+
  xlab("Mean and Median branch length values") +
  ylab("Density")
  
ggsave("Plot_for_non_outliers.pdf", plot = plot,width = 7 , height = 5 ,dpi = 1200)