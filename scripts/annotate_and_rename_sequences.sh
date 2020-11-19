#!/bin/bash

#This script  as an argument takes fna file folder and location of rename_fasta.py script
#Example annotate_and_rename_sequences.sh Streptomyces.csv

FOLDER="$1"
FILENAME="$(basename -s .csv $FOLDER)"
FOLDER_2="$2"
MSA_ALG=$3
TREE_ALG=$4
SPLIT=$5

hr() {
  local start=$'\e(0' end=$'\e(B' line='qqqqqqqqqqqqqqqq'
  local cols=${COLUMNS:-$(tput cols)}
  while ((${#line} < cols)); do line+="$line"; done
  printf '%s%s%s\n' "$start" "${line:0:cols}" "$end"
}

if [ "$SPLIT" -eq 0 ]; then
	echo "Step 1: Downloading genomes..."
	cat $FOLDER | awk -F ',' '{print $15}' | sed 's/ftp:/rsync:/g' | sed 's/"//g' | tail -n +2  > list_tmp.txt

	cat list_tmp.txt | awk -F '/' '{print $10}' |  sed 's/"//g' |sed 's/$/_genomic.fna.gz/'  > suffix.txt

	paste -d '/' list_tmp.txt suffix.txt > to_download.txt

	rm -f list_tmp.txt suffix.txt

	cat to_download.txt |parallel -j 25 rsync --copy-links --times --verbose {} ./

	ls *.fna.gz | parallel gunzip {}

	rm to_download.txt
	hr

	echo "Step 2: Running barrnap..."
	ls *.fna | parallel barrnap --quiet --outseq {.}.rna {}
	hr
	echo "Step 3: Extracting 16S rRNAs..."
	for i in *.rna; do  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | grep -A 1 "^>16S" > $i.16S ; done
	hr
else
	echo "Step 0: Splitting genomes..."
	split -l $SPLIT -a 2 -d $FOLDER
	for i in x*;do { echo '#Organism Name,Organism Groups,Strain,BioSample,BioProject,Assembly,Level,Size(Mb),GC%,Replicons,WGS,Scaffolds,CDS,Release Date,GenBank FTP,RefSeq FTP'; cat $i; } > $i.csv ;done
	rm x00.csv
	mv x00 x00.csv
	for i in x*.csv; do
		echo "Step 1: Downloading genomes..."
		cat $i | awk -F ',' '{print $15}' | sed 's/ftp:/rsync:/g' | sed 's/"//g' | tail -n +2  > list_tmp.txt

		cat list_tmp.txt | awk -F '/' '{print $10}' |  sed 's/"//g' |sed 's/$/_genomic.fna.gz/'  > suffix.txt

		paste -d '/' list_tmp.txt suffix.txt > to_download.txt

		rm -f list_tmp.txt suffix.txt

		cat to_download.txt |parallel -j 25 rsync --copy-links --times --verbose {} ./

		ls *.fna.gz | parallel gunzip {}

		rm to_download.txt
		hr

		echo "Step 2: Running barrnap..."
		ls *.fna | parallel barrnap --quiet --outseq {.}.rna {}
		hr
		echo "Step 3: Extracting 16S rRNAs..."
		for i in *.rna; do  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i | grep -A 1 "^>16S" > $i.16S ; done
		hr
		rm -f *.fna
	done
fi
echo "Step 4: Renaming sequence files..."
#pass arguments as follow: rename_fasta.py /datasets/Streptomyces.csv 
python "${FOLDER_2}"rename_fasta.py $FOLDER
mkdir "${FILENAME}_no_rrna"
mv renamed_rna $FILENAME

find -size 0 -print > list_2.txt
while read p; do cp $p "${FILENAME}_no_rrna" ; done < list_2.txt

rm list_2.txt


cd "${FILENAME}_no_rrna" && rm -f *.16S list_2.txt && cd ..

cd $FILENAME
hr
echo "Step 5: Removing 1 sequence files and zero size files..."
mkdir one_rrna
mv `grep -c  '>' *.fasta | awk -F':' '{ if ($2<=1) print $1}'` one_rrna

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
	ls *.fasta.mafft | parallel 'iqtree --quiet  -T 1 -m GTR -s {}'
	for i in *.treefile; do cat $i > $i.nwk; done
elif [ "$TREE_ALG" = raxml ]; then
	ls *.fasta.mafft | parallel 'raxmlHPC -s {} -n {.}.tmp -m  GTRCAT --print-identical-sequences'
	for i in RAxML_bestTree.*.tmp; do cat $i > $i.nwk; done
else
	echo "Named TREE algorithm was not found! "
	echo "Please check the spelling and look in help for available options"
	exit 1
fi

mv *.nwk ../ 
mv *.fasta ../ 
mv one_rrna ../
cd ../

rm `find ./ -size 0 | grep "nwk"`
hr
echo "Step 8: Computing dataframes out of branch lengths..."
Rscript --vanilla "${FOLDER_2}"compute_dataframes.R $FILENAME ./

rm -f *.nwk *.pdf *.fna *.fai *.rna *.16S
rm -rf $FILENAME

mkdir "${FILENAME}_16S_rna"
cd one_rrna && basename -s .fasta `ls *` > "${FILENAME}_one_tmp.csv" && { echo 'Names'; cat "${FILENAME}_one_tmp.csv"; } > "${FILENAME}_one.csv" && mv "${FILENAME}_one.csv" ../ && rm *.csv && cd ..
mv one_rrna "${FILENAME}_one_rna" 
mv *.fasta  "${FILENAME}_16S_rna"