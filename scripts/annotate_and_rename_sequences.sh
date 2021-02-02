#!/bin/bash

# Passed arguments:
# $1 - The name of a database to build
# $2 - Th folder, where are all teh scripts located
# $3 - MSA algorithm
# $4 - Phylogeny algorithm
# $5 - Value of how many genomes split download

# Set script variabels
FOLDER="$1"
FILENAME="$(basename -s .csv $FOLDER)"
FOLDER_2="$2"
MSA_ALG=$3
TREE_ALG=$4
SPLIT=$5

# Set drawing line function
hr() {
  local start=$'\e(0' end=$'\e(B' line='qqqqqqqqqqqqqqqq'
  local cols=${COLUMNS:-$(tput cols)}
  while ((${#line} < cols)); do line+="$line"; done
  printf '%s%s%s\n' "$start" "${line:0:cols}" "$end"
}

# Split logic for genome download
if [ "$SPLIT" -eq 0 ]; then
	# No split. Download in one big chunk
	echo "Step 1: Downloading genomes..."
	# List all the rsync:// links to use
	cat $FOLDER | awk -F ',' '{print $15}' | sed 's/ftp:/rsync:/g' \
	| sed 's/"//g' | tail -n +2  > list_tmp.txt
	# Create suffix list (to download only genome sequences) add there _genomic.fna.gz
	cat list_tmp.txt | awk -F '/' '{print $10}' |  sed 's/"//g' \
	|sed 's/$/_genomic.fna.gz/'  > suffix.txt

	# Combine the rsync:// and suffix lists to final list with download links
	paste -d '/' list_tmp.txt suffix.txt > to_download.txt

	# Remove previous lists
	rm -f list_tmp.txt suffix.txt

	# Downloading the files....
	cat to_download.txt |parallel -j 25 rsync --copy-links --times --verbose {} ./

	# Extract them from archives
	ls *.fna.gz | parallel gunzip {}

	# Final cleaning
	rm to_download.txt
	hr

	# Running barrnap for 16S rRNA annotation
	echo "Step 2: Running barrnap..."
	ls *.fna | parallel barrnap --quiet --outseq {.}.rna {}
	hr
	# Extract only 16S sequences out of annotated files
	echo "Step 3: Extracting 16S rRNAs..."
	for i in *.rna; do \
	 awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i \
	 | grep -A 1 "^>16S" > $i.16S ; done
	hr
else
	# TO-DO: first do the cleaning of the csv, then split the link list. 
	# In that way you should not care about the missing header


	# Need to split the download (So we would download, then extract from archives
	# then annonate with barrnap, then delete the genomes and repeat the download)
	echo "Step 0: Splitting genomes..."
	# Split the .csv file based on number of rows (organisms to download)
	split -l $SPLIT -a 2 -d $FOLDER
	# Add the header to newly generated (from split) csv files
	for i in x*;do { echo '#Organism Name,Organism Groups,Strain,BioSample,BioProject,Assembly,Level,Size(Mb),GC%,Replicons,WGS,Scaffolds,CDS,Release Date,GenBank FTP,RefSeq FTP'; cat $i; } > $i.csv ;done
	# File x00 now contains 2 headers. So we remove the modified one, and use the 
	# one from splitting
	rm x00.csv
	mv x00 x00.csv

	#Download, extract, annotate, delete process for every csv file
	for i in x*.csv; do
		echo "Step 1: Downloading genomes..."

		# The below steps are the same as described above, so:
		# Create clean link file for downloading
		cat $i | awk -F ',' '{print $15}' | sed 's/ftp:/rsync:/g' | sed 's/"//g' \
		| tail -n +2  > list_tmp.txt

		cat list_tmp.txt | awk -F '/' '{print $10}' |  sed 's/"//g' \
		|sed 's/$/_genomic.fna.gz/'  > suffix.txt

		paste -d '/' list_tmp.txt suffix.txt > to_download.txt
		rm -f list_tmp.txt suffix.txt
		# Download and extract files
		cat to_download.txt |parallel -j 25 rsync --copy-links --times --verbose {} ./
		ls *.fna.gz | parallel gunzip {}
		rm to_download.txt
		hr

		# Run barrnap on those files 
		echo "Step 2: Running barrnap..."
		ls *.fna | parallel barrnap --quiet --outseq {.}.rna {}
		hr
		# Extract 16S rRNAs from those files
		echo "Step 3: Extracting 16S rRNAs..."
		for i in *.rna; do \
		 awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i \
		 | \grep -A 1 "^>16S" > $i.16S ; done
		hr
		rm -f *.fna
	done
fi

# Rename all the 16S files, using csv filewe used for download.
echo "Step 4: Renaming sequence files..."
python "${FOLDER_2}"rename_fasta.py $FOLDER
mkdir "${FILENAME}_no_rrna"
mv renamed_rna $FILENAME

# Find size 0 files (no rrnas or 16S rrnas). Move them to the separate directory
find -size 0 -print > list_2.txt
while read p; do cp $p "${FILENAME}_no_rrna" ; done < list_2.txt
rm list_2.txt
cd "${FILENAME}_no_rrna" && rm -f *.16S list_2.txt && cd ..
# Go to the directory with renamed files
cd $FILENAME
hr

# Moving one sequence files to a separate directory. They will not be included
# in further analysis. They are used in the phylogeny later with --tree argument
echo "Step 5: Removing 1 sequence files and zero size files..."
mkdir one_rrna
mv `grep -c  '>' *.fasta | awk -F':' '{ if ($2<=1) print $1}'` one_rrna

rm `find ./ -size 0 | grep "fasta"`


# Add number to the duplicate fasta headers in fasta files
for i in *.fasta; do cat $i | perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; '\
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

# Remove temporary files
mv *.nwk ../ 
mv *.fasta ../ 
mv one_rrna ../
cd ../

# Remove zero bytes .nwk files. Actually they points to the failed phylogeny run
rm `find ./ -size 0 | grep "nwk"`
hr

# Compute dataframe statistics
echo "Step 8: Computing dataframes out of branch lengths..."
Rscript --vanilla "${FOLDER_2}"compute_dataframes.R $FILENAME ./

# Remove all the unnesecary files
rm -f *.nwk *.pdf *.fna *.fai *.rna *.16S
rm -rf $FILENAME

mkdir "${FILENAME}_16S_rna"

# Make a dataframe one 1 seq fasta sequences
cd one_rrna && basename -s .fasta `ls *` > "${FILENAME}_one_tmp.csv" && \
{ echo 'Names'; cat "${FILENAME}_one_tmp.csv"; } > "${FILENAME}_one.csv" \
&& mv "${FILENAME}_one.csv" ../ && rm *.csv && cd ..

# Final cleaning 
mv one_rrna "${FILENAME}_one_rrna" 
mv *.fasta  "${FILENAME}_16S_rna"