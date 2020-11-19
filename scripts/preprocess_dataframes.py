# Passed arguments:
# $1 - The input dataframe file. Should be in a format as precomputed database 
#      one
# $2 - The folder where write the results
# $3 - Directory from where execute this script

import pandas as pd
import sys
import os

# Set working directory
os.chdir(sys.argv[3])

# Read the dataframe, where non outliers are keeped. Filters the fasta to remain 
# only one 16S sequence
to_filter = pd.read_csv(str(sys.argv[1])+"_no_ouliers.csv")
# Make a list out of Names to iterate through
to_iterate = list(to_filter['Names'])

# Generate remain_one.txt file, where organims, which are non outliers are 
# listed.
with open(sys.argv[2]+"/remain_one.txt", "w") as f:
  f.write('\n'.join(to_iterate))
  f.write('\n')

# Read the mean outlier dataframe and create txt file with the organisms
to_filter = pd.read_csv(str(sys.argv[1])+"_mean_outliers.csv")
to_iterate = list(to_filter['Names'])
with open(sys.argv[2]+"/outliers_mean.txt", "w") as f:
  f.write('\n'.join(to_iterate))
  f.write('\n')

# Read th median outlier dataframe and generate a txt file where they are listed
to_filter = pd.read_csv(str(sys.argv[1])+"_median_outliers.csv")
to_iterate = list(to_filter['Names'])
with open(sys.argv[2]+"/outliers_median.txt", "w") as f:
  f.write('\n'.join(to_iterate))
  f.write('\n')

# Read the dataframe with organisms who have only one 16S rRNA. List them in txt
# file
to_filter = pd.read_csv(str(sys.argv[1])+"_one.csv")
to_iterate = list(to_filter['Names'])
with open(sys.argv[2]+"/only_one.txt", "w") as f:
  f.write('\n'.join(to_iterate))
  f.write('\n')