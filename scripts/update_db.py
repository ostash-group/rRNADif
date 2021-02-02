# $1 - database name (database to update)
# $2 - file to use for update
# $3 - treat $2 as file or as folder with sequences? ("TRUE" = folder)
# $4 - working directory


import pandas as pd
import os
import sys

os.chdir(sys.argv[4])


def process_dataframe(file):
    '''
    Returns processed dataframe from standard csv file, as downloaded from NCBI
    The dataframe contains two colums Code and Organism
    Code is a filename, as the last /string in GenBank FTP field
    Organism is an organism name, as name + strains, if strain is not listed in name.
    '''
    #Need to properly display all names. If less - strings are incomplete
    pd.options.display.max_colwidth = 100
    df = pd.read_csv(file)
    #Make tmp dataframe with splited GenBank FNT field
    df = df[df["GenBank FTP"].notna()]
    df['id'] = list(range(len(df)))
    df1 = df["GenBank FTP"].str.split(pat = '/')
    #Make new dataframe to append to
    df_clean = pd.DataFrame(columns = ['Code', 'Organism', 'id'])
    for i,row in enumerate(df1.iteritems()):
        #Checks if information in strain field is already in name field,
        #If not, concatenate name and atrin fileds
        if df['#Organism Name'].iloc[i].split(' ')[-1] == str(str(df['Strain'].iloc[i]).split(' ')[-1]):
            name = df['#Organism Name'].iloc[i]
            name = name.replace(' ', '_').replace('/', '_').replace(':', '_')\
          .replace(';', '_').replace(',', '_').replace('[','').replace(']','') 
        else:
            name = df['#Organism Name'].iloc[i] +' '+str(df['Strain'].iloc[i]).replace("/",'_')
            name = name.replace(' ', '_').replace('/', '_').replace(':', '_')\
          .replace(';', '_').replace(',', '_').replace('[','').replace(']','') 
        #Appends information to clean dataframe
        df_clean =  df_clean.append(pd.Series([row[-1][-1], name, i], index = df_clean.columns), 
        ignore_index = True)
    df_clean_1 = pd.merge(df, df_clean, on='id')
    return df_clean_1

def update_dataframe(master, subset):
    '''
    Update master dataframe, using info from the subset one (which entries 
    are not present there)
    The master dataframe can be the one with 'Species_names' column (then it
    will be used), or from "Browse by organism" csv file. (then preprocessing
    will be done)
    '''
    
    if 'Species_names' in master.columns:
      subset_final_1 = subset[~subset.Names.isin(list(master['Species_names']))]
      subset_final = subset_final_1['Filename']
    else:
      subset_final = subset[~subset.Names.isin(list(master['Organism']))]

    return subset_final

# Main program
# Read the data
master_csv = pd.read_csv("../datasets/downloaded_csv/"+sys.argv[1]+".csv")
# If it is a folder data -> preprocess
if sys.argv[3] == "TRUE":
  subset_csv = pd.DataFrame(columns = ['Names', 'id', 'Filename']) 
  with open("list.txt") as f:
    files = f.readlines()
    for i, f_l in enumerate(files):
      fl, fe = os.path.splitext(f_l)
      subset_csv = subset_csv.append(pd.Series([fl, i, f_l.rstrip()], index = subset_csv.columns), ignore_index=True)
  
else:
# Else use the function to preprocess dataframe
  subset_csv = process_dataframe(sys.argv[2])
  subset_csv.rename({'Organism':"Names"},axis=1, inplace=True)

to_write = update_dataframe(master_csv,subset_csv)
# Write results
to_write.to_csv("to_update.csv", index = False)