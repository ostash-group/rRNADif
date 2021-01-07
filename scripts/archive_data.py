# $1 - dataset to transorm. Must be "Browse by Genome" csv file
# $2 - working directory

import pandas as pd
import sys
import os

os.chdir(sys.argv[2])

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
    df_clean = pd.DataFrame(columns = ['Organism'])
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
        df_clean =  df_clean.append(pd.Series( [name], index = df_clean.columns), 
        ignore_index = True)
    return df_clean

# Process the dataframe. Store the result
data = process_dataframe(sys.argv[1])
filename, trash = os.path.splitext(str(sys.argv[1]))
data.to_csv(str(filename) + "_new.csv", index = False)