# Passed arguments:
# $1 - filename for a file to open 
# $2 - folder, where to write the newly generated file (with single sequence
# $3 - folder, where to find initial file, to match the $1 argument 

from Bio import SeqIO
import sys

# Sequqnce counter
i = 1
# Open file. Save the id and seq in the variables. Increment i, if >1 -> exit
for seq in SeqIO.parse(sys.argv[3]+"/"+sys.argv[1] + ".fasta", "fasta"):
  if i == 1:
    ID = seq.id
    SEQ = seq.seq
    i+=1
  else:
    exit

# Write this ID and SEQ variables into file
with open(sys.argv[2]+ "/" +str(ID)+"_one_seq"+".clean", "w") as f:
  f.write('>'+str(ID)+'\n'+str(SEQ) + '\n')