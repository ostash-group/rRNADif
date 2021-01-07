#!/bin/bash

# Passed arguments:
# $1 - The filename of a folder
# $2 - Where to locate R and python scripts
# $3 - Which MSA algoritm to use
# $4 - Which phylogeny computign program to use

# Set the arguments to the script variables
FILENAME="$1"
FOLDER_2="$2"
MSA_ALG=$3
TREE_ALG=$4

# Define the drawing line function
hr() {
  local start=$'\e(0' end=$'\e(B' line='qqqqqqqqqqqqqqqq'
  local cols=${COLUMNS:-$(tput cols)}
  while ((${#line} < cols)); do line+="$line"; done
  printf '%s%s%s\n' "$start" "${line:0:cols}" "$end"
}

# Run the barrnap
echo "Step 2: Running barrnap..."
cat list.txt | parallel barrnap --quiet --outseq {.}.rna {}
hr

# Extract 16S rRNA. Here we have some grep magic. 
echo "Step 3: Extracting 16S rRNAs..."
for i in *.rna; do  \
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i \
| grep -A 1 "^>16S" > $(basename -s .rna $i).16S ; done
hr

# Rename the sequence headers and files, using filename.
echo "Step 4: Renaming sequence files..."
for i in *.16S ; do python ../scripts/rename_fasta_headers.py $(pwd) $i
done

for i in *.16S ; do 
perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' $i > $(basename -s .16S $i).fasta
done

# Find size 0 files (no rrnas or 16S rrnas). Move them to the separate directory
mkdir "${FILENAME}_no_rrna"
mkdir "${FILENAME}"
mv *.fasta $FILENAME
find -size 0 -print > list_2.txt
while read p; do cp $p "${FILENAME}_no_rrna" ; done < list_2.txt

# Remove some temporary files
rm list_2.txt
rm "${FILENAME}_no_rrna"/list_2.txt
cd $FILENAME
hr

# Moving one sequence files to a separate directory. They will not be included
# in further analysis. They are used in the phylogeny later with --tree argument
echo "Step 5: Removing 1 sequence files and zero size files..."
mkdir one_rrna
mv `grep -c  '>' *.fasta | awk -F':' '{ if ($2<=1) print $1}'` one_rrna

#Remove size 0 fasta files 
rm `find ./ -size 0 | grep "fasta"`

# Add number to the duplicate fasta headers in fasta files
for i in *.fasta; do cat $i | perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' \
| sed 's/;_/_/g' > $i.clean; done
hr

# Perform the MSA with the chosen algorithm 
echo "Step 6: Performing mafft alignment..."
if [ "$MSA_ALG" = clustalo ]; then
	ls *.clean | parallel  'clustalo -i {} -o {.}.mafft'
elif [ "$MSA_ALG" = muscle ]; then
	ls *.clean | parallel 'muscle -quiet  -in {} -out {.}.mafft'
elif [ "$MSA_ALG" = mafft ]; then
	ls *.clean | parallel  'mafft --quiet {} > {.}.mafft'
else
	echo "Named MSA algorithm was not found! "
	echo "Please check the spelling and look in help for available options"
	exit 1
fi
hr

# Perform the phylogeny computation with the chosen algoritm
echo "Step 7: Building phylogenetic trees..."
if  [ "$TREE_ALG" = fasttree ]; then
	ls *.fasta.mafft | parallel 'fasttree -quiet -nt -gtr -out {.}.nwk  {}'
elif [ "$TREE_ALG" = iqtree ]; then
	ls *.fasta.mafft | parallel 'iqtree --quiet -T 1 -m GTR -s {}'
	for i in *.treefile; do cat $i > $i.nwk; done
elif [ "$TREE_ALG" = raxml ]; then
	ls *.fasta.mafft | parallel 'raxmlHPC -s {} -n {.}.tmp -m  GTRCAT --print-identical-sequences'
	for i in RAxML_bestTree.*.tmp; do cat $i > $i.nwk; done
else
	echo "Named TREE algorithm was not found! "
	echo "Please check the spelling and look in help for available options"
	exit 1
fi
ls *.fasta.mafft | parallel 'fasttree -quiet -nt -gtr -out {.}.nwk  {}'

# Remove all unnecessary files
mv *.nwk ../ 
mv *.fasta ../ 
mv one_rrna ../
cd ../

# Remove size 0 phylogenetic files (files, where previous step failes)
rm `find ./ -size 0 | grep "nwk"`
hr

# Compute dataframe statistics
echo "Step 8: Computing dataframes out of branch lengths..."
Rscript --vanilla ../scripts/compute_dataframes.R $FILENAME $(pwd)

# Remove all the unnesecary files
rm -f *.nwk *.pdf *.fna *.fai *.rna *.16S *.txt *.sh 
rm -rf $FILENAME

mkdir "${FILENAME}_16S_rna"

# Make a dataframe one 1 seq fasta sequences
cd one_rrna && basename -s .fasta `ls *` > "${FILENAME}_one_tmp.csv" && \
{ echo 'Names'; cat "${FILENAME}_one_tmp.csv"; } > "${FILENAME}_one.csv" && \
mv "${FILENAME}_one.csv" ../ && rm *.csv && cd ..

# Final cleaning 
mv one_rrna "${FILENAME}_one_rrna" 
mv *.fasta  "${FILENAME}_16S_rna"
