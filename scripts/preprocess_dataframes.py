# Passed arguments:
# $1 - The input dataframe file. Should be in a format as precomputed database 
#      one
# $2 - The folder where write the results
# $3 - Directory from where execute this script
# $4 - 0 if use all strains, 1 - if use only the first strains of a species 
#      listed

import pandas as pd
import sys
import os



used_names = []


def only_species(to_iterate):
    global used_names
    to_iterate_new = []
    for name in to_iterate:
        name_list = name.split('_')
        new_name = name_list[0]+'_'+name_list[1]
        if new_name in used_names:
          continue
        else:
          used_names.append(new_name)
          to_iterate_new.append(name)
    return to_iterate_new

# Set working directory
os.chdir(sys.argv[3])

# Read the dataframe, where non outliers are keeped. Filters the fasta to remain 
# only one 16S sequence
to_filter = pd.read_csv(str(sys.argv[1])+"_no_ouliers.csv")
# Make a list out of Names to iterate through
to_iterate = list(to_filter['Names'])

# Generate remain_one.txt file, where organims, which are non outliers are 
# listed.
if sys.argv[4] == "TRUE":
	to_iterate_new = only_species(to_iterate)
else:
	to_iterate_new = to_iterate
if len(to_iterate_new) != 0:
  with open(sys.argv[2]+"/remain_one.txt", "w") as f:
    f.write('\n'.join(to_iterate_new))
    f.write('\n')

# Read the mean outlier dataframe and create txt file with the organisms
to_filter = pd.read_csv(str(sys.argv[1])+"_mean_outliers.csv")
to_iterate = list(to_filter['Names'])
if sys.argv[4] == "TRUE":
	to_iterate_new = only_species(to_iterate)
else:
	to_iterate_new = to_iterate

if len(to_iterate_new) != 0:
  with open(sys.argv[2]+"/outliers_mean.txt", "w") as f:
    f.write('\n'.join(to_iterate_new))
    f.write('\n')

# Read th median outlier dataframe and generate a txt file where they are listed
to_filter = pd.read_csv(str(sys.argv[1])+"_median_outliers.csv")
to_iterate = list(to_filter['Names'])
if sys.argv[4] == "TRUE":
	to_iterate_new = only_species(to_iterate)
else:
	to_iterate_new = to_iterate
if len(to_iterate_new) != 0:
  with open(sys.argv[2]+"/outliers_median.txt", "w") as f:
    f.write('\n'.join(to_iterate_new))
    f.write('\n')

# Read the dataframe with organisms who have only one 16S rRNA. List them in txt
# file
to_filter = pd.read_csv(str(sys.argv[1])+"_one.csv")
to_iterate = list(to_filter['Names'])
if sys.argv[4] == "TRUE":
	to_iterate_new = only_species(to_iterate)
else:
	to_iterate_new = to_iterate
if len(to_iterate_new) != 0:
  with open(sys.argv[2]+"/only_one.txt", "w") as f:
    f.write('\n'.join(to_iterate_new))
    f.write('\n')
