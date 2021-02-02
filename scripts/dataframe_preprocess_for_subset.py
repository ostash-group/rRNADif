# Passed arguments:
# $1 - The input database file. Used for subsetting the master precomputed one
# $2 - A file (path) which to write as a subset of master database
# $3 - The folder where to look for a precomputed database
# $4 - The database name to gather the _all and _one files
# $5 - A file (path) which to write as a subset of "one" database. So we gather
#      the organisms that have only one 16S in their genome 


# Import linraries
import pandas as pd
import sys
import os

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
    df1 = df["GenBank FTP"].str.split(pat = '/')
    #Make new dataframe to append to
    df_clean = pd.DataFrame(columns = ['Code', 'Organism'])
    i = 0
    for row in df1.iteritems():
        #Checks if information in strain field is already in name field,
        #If not, concatenate name and atrin fileds
        if df['#Organism Name'].iloc[i].split(' ')[-1] == str(str(df['Strain'].iloc[i]).split(' ')[-1]):
            name = df['#Organism Name'].iloc[i]
        else:
            name = df['#Organism Name'].iloc[i] +' '+str(df['Strain'].iloc[i]).replace("/",'_')
        #Appends information to clean dataframe
        df_clean =  df_clean.append(pd.Series([row[-1][-1], name], index = df_clean.columns), 
        ignore_index = True)
        i+=1
    return df_clean


def check_two_dataframes(master, subset):
    '''
    Compare to dataframes - master one (compare to) and subset one (which to 
    compare).
    The master dataframe should contain 'Organism' column. The subset is a 
    csv from Genbank, organisms from which are going to be renamed and checked 
    with the master 'Organism' column
    '''
    cleaned_subset =  [name.replace(' ', '_').replace('/', '_').replace(':', '_')\
        .replace(';', '_').replace(',', '_').replace('[','').replace(']','') \
            for name in list(subset['Organism'])]
    subset_final = master[master.Names.isin(cleaned_subset)]
    return subset_final

# Read master dataframe files and subset one
master_csv = pd.read_csv(sys.argv[3]+"/" +sys.argv[4]+"_all.csv")
master_2_csv = pd.read_csv(sys.argv[3]+"/" +sys.argv[4]+"_one.csv")
subset_csv = process_dataframe(sys.argv[1])
# Use the dataframe with organisms that have > 2 16S for a subset
subset_final = check_two_dataframes(master_csv, subset_csv)
subset_final.to_csv(sys.argv[2], index = False)
# Use the dataframe with organisms that have one 16S for a subset
subset_final_2 = check_two_dataframes(master_2_csv, subset_csv)
subset_final_2.to_csv(sys.argv[5], index = False)