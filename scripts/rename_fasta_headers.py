# Script renames fasta headers, using filename
# $1 - working directory
# $2 - file to rename
import sys
import os
from Bio import SeqIO

# Change the working directory
os.chdir(sys.argv[1])

# Open the files, get the filename and write them
with open(sys.argv[2], "r") as f:
    seqs = []
    filename, f_ext = os.path.splitext(sys.argv[2])
    for record in SeqIO.parse(f, "fasta"):
        seqs.append(record.seq)

with open(sys.argv[2], "w") as f_new:
    for i in range(len(seqs)):
        f_new.write(">"+ filename +"\n"+ str(seqs[i])+"\n")