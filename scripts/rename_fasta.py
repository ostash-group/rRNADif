import pandas
import numpy as np
import sys
import os
import re
from Bio import SeqIO
#Globals
seqs = []
names = []
already_seen = []


def set_path(string, new_folder_string):
    '''
    Sets path to desired place.
    Returns two variables path and new_path (with folder string)
    New path can be used to store result files
    NEW FOLDER STRING MUST BE WITH SLASH
    '''
    os.chdir(string)
    path = os.getcwd()
    new_path = path + new_folder_string
    #Creates new folder with given string
    try:
        os.mkdir(new_path)
    except:
        print("New folder was not created. Mybe it already exists")
    return path, new_path


def get_files_in_dir(path=os.getcwd(), extension = False):
    '''
    Returns a list with listed files from a path with given extension
    If extension is not given, listing all files
    Returns a list
    '''
    if extension:
        entries = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) and re.search(r""+extension, os.path.splitext(f)[1])]
    else:
        entries = [f for f in os.listdir(path) ]
    return entries


def get_organism_info(filename, dataframe):
    '''
    Get organism name from processed csv file using filename as a code
    '''
    #Returns pandas series only Organism field, where Code field matches filename, converts it to string,
    #and cut leading zeroes
    return dataframe[dataframe['Code']==filename]['Organism'].to_string(header = False, index = False).lstrip()


def get_organism_code(name, dataframe):
    '''
    Returns organism code, given it's name as a code
    '''
    #Returns pandas series only Code field, where Organism field matches filename, converts it to string,
    #and cut leading zeroes
    return dataframe[dataframe['Organism']==name]['Code'].to_string(header = False, index = False)


def process_dataframe(file):
    '''
    Returns processed dataframe from standard csv file, as downloaded from NCBI
    The dataframe contains two colums Code and Organism
    Code is a filename, as the last /string in GenBank FTP field
    Organism is an organism name, as name + strains, if strain is not listed in name.
    '''
    #Need to properly display all names. If less - strings are incomplete
    pandas.options.display.max_colwidth = 100
    df = pandas.read_csv(file)
    #Make tmp dataframe with splited GenBank FNT field
    df = df[df["GenBank FTP"].notna()]
    df1 = df["GenBank FTP"].str.split(pat = '/')
    #Make new dataframe to append to
    df_clean = pandas.DataFrame(columns = ['Code', 'Organism'])
    i = 0
    for row in df1.iteritems():
        #Checks if information in strain field is already in name field,
        #If not, concatenate name and atrin fileds
        if df['#Organism Name'].iloc[i].split(' ')[-1] == str(str(df['Strain'].iloc[i]).split(' ')[-1]):
            name = df['#Organism Name'].iloc[i]
        else:
            name = df['#Organism Name'].iloc[i] +' '+str(df['Strain'].iloc[i]).replace("/",'_')
        #Appends information to clean dataframe
        df_clean =  df_clean.append(pandas.Series([row[-1][-1], name], index = df_clean.columns), ignore_index = True)

        i+=1
    return df_clean


def rename_fasta_headers(record):
    '''
    Checks if given record contain "16S" in it's description
    If yes, appends sequence to the global seqs list
    '''
    global seqs
    if len(re.findall(r"16S", record.description)):
        seqs.append(str(record.seq))


def add_fasta_headers(record):
    '''
    Just adds entries from fasta file to list.
    Need for proper work of write.fasta() function. List then used as input there
    '''
    global seqs
    seqs.append(str(record.seq))


def write_fasta():
    '''
    Writes fasta file with organism names and seqs. Returns clean seqs and names list for next entry
    If there is no 16S listed in the file, add filename to add_rnas and organism name to add_rnas_names
    '''
    global seqs, names, new_path, df, add_rnas,add_rnas_names
    if check_if_are_entries():
        name_tmp = names[0].replace(' ', '_').replace('/', '_').replace(':', '_').replace(';', '_').replace(',', '_').replace('[','').replace(']','').replace('*','').replace("'", "")
        if name_tmp not in already_seen:
            with open(new_path +'/' +name_tmp+".fasta", "w") as f:
                for i in range(len(seqs)):
                    f.write(">" + name_tmp + "\n" +seqs[i] + "\n")
            already_seen.append(name_tmp)
        else:
            for i in range(1,10000):
                if name_tmp+"__"+str(i) in already_seen:
                    pass
                else:
                    name_tmp_1 = name_tmp + "__" + str(i)
                    with open(new_path +'/' +name_tmp_1+".fasta", "w") as f:
                        for i in range(len(seqs)):
                            f.write(">" + name_tmp_1 + "\n" +seqs[i] + "\n")
                    already_seen.append(name_tmp_1)
                    break
    else:
        name_tmp = names[0].replace(' ', '_').replace('/', '_').replace(':', '_').replace(';', '_').replace(',', '_').replace('[','').replace(']','').replace('*','').replace("'", "")
        with open(new_path +'/' +name_tmp+".fasta", "w") as f:
            f.write("")

def zero_state():
    '''
    Returns variables to it's initial state
    '''
    global seqs, names
    seqs = []
    names = []


def check_if_are_entries():
    '''
    Checks if there is any entries in given list
    And we are using this list to append found 16S entries. 
    So if it is empty - need to recheck manually
    Returns True if there are sequences, False if none
    '''
    global seqs, add_rnas, add_rnas_names
    if len(seqs) == 0 :
        add_rnas.append(get_organism_code(names[0],df))
        add_rnas_names.append(names[0])
        return False
    else:
        return True


path, new_path = set_path(os.getcwd(),'/renamed_rna')
df = process_dataframe(sys.argv[1])
entries = get_files_in_dir(path, "16S")

print("Renaming " + str(len(entries)) + " files in a directory")

add_rnas = []
add_rnas_names = []
for entry in entries:
    file = open(entry, 'r')
    filename, file_extension = os.path.splitext(entry)
    name = get_organism_info(filename.replace('_genomic.rna',''), df)
    names.append(name.lstrip())
    for record in SeqIO.parse(file, 'fasta'):
        rename_fasta_headers(record)
    write_fasta()
    zero_state()