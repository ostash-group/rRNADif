#!/bin/bash

#This script  as an argument takes fna file folder and location of rename_fasta.py script
#Example annotate_and_rename_sequences.sh Streptomyces.csv

FILENAME="$1"
FOLDER_2="$2"
MSA_ALG=$3
TREE_ALG=$4

hr() {
  local start=$'\e(0' end=$'\e(B' line='qqqqqqqqqqqqqqqq'
  local cols=${COLUMNS:-$(tput cols)}
  while ((${#line} < cols)); do line+="$line"; done
  printf '%s%s%s\n' "$start" "${line:0:cols}" "$end"
}


echo "Step 2: Running barrnap..."
cat list.txt | parallel barrnap --quiet --outseq {.}.rna {}
hr
echo "Step 3: Extracting 16S rRNAs..."
for i in *.rna; do  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | grep -A 1 "^>16S" > $(basename -s .rna $i).16S ; done
hr
echo "Step 4: Renaming sequence files..."
#pass arguments as follow: rename_fasta.py /datasets/Streptomyces.csv 
for i in *.16S; do awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-6); next} 1' $i | perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' > $(basename -s .16S $i).fasta; done
mkdir "${FILENAME}_no_rrna"
mkdir "${FILENAME}"
mv *.fasta $FILENAME
find -size 0 -print > list_2.txt

while read p; do cp $p "${FILENAME}_no_rrna" ; done < list_2.txt

rm list_2.txt

rm "${FILENAME}_no_rrna"/list_2.txt

cd $FILENAME
hr
echo "Step 5: Removing 1 sequence files and zero size files..."
mkdir one_rrna
mv `grep -c  '>' *.fasta | awk -F':' '{ if ($2<=1) print $1}'` one_rrna

rm `find ./ -size 0 | grep "fasta"`

for i in *.fasta; do cat $i | perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' | sed 's/;_/_/g' > $i.clean; done
hr
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


mv *.nwk ../ 
mv *.fasta ../ 
mv one_rrna ../
cd ../

rm `find ./ -size 0 | grep "nwk"`
hr
echo "Step 8: Computing dataframes out of branch lengths..."

Rscript --vanilla ../scripts/compute_dataframes.R $FILENAME $(pwd)
rm -f *.nwk *.pdf *.fna *.fai *.rna *.16S *.txt *.sh 
rm -rf $FILENAME

mkdir "${FILENAME}_16S_rna"
cd one_rrna && basename -s .fasta `ls *` > "${FILENAME}_one_tmp.csv" && { echo 'Names'; cat "${FILENAME}_one_tmp.csv"; } > "${FILENAME}_one.csv" && mv "${FILENAME}_one.csv" ../ && rm *.csv && cd ..
mv one_rrna "${FILENAME}_one_rna" 
mv *.fasta  "${FILENAME}_16S_rna"