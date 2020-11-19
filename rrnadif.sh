#!/bin/bash
#
# AUTHOR: Pavlo Hrab

# Global variables
VERSION=0.1
# The default is precomputed, but the -n flag rewrites this. Points to the 
# folder
DATABASE_FOLDER=Bacteria_full
PLOT=0
TREE=0
INPUT_BOOL=0
NAME_BOOL=0
DATABASE_BOOL=0
MAKEDB_BOOL=0
MAKEDB_BOOL_1=0
STEP=0
REDO=0
SEQUENCE_B=0
MSA_ALG=mafft
TREE_ALG=fasttree
PHYLO=0
SPLIT=100
SPLIT_BOOL=0

# If any errors stop script execution
set -e 

# Define the divider function. It's dynamically goes to the current
# width of the terminal
hr() {
  local start=$'\e(0' end=$'\e(B' line='qqqqqqqqqqqqqqqq'
  local cols=${COLUMNS:-$(tput cols)}
  while ((${#line} < cols)); do line+="$line"; done
  printf '%s%s%s\n' "$start" "${line:0:cols}" "$end"
}

hr

# Define the help message
help(){
cat << EOF

rRNADif  Version: ${VERSION}

rRNADif computes intragenomic 16S rRNA variability for provided genome sequence 
against precomputed database, using branch lenghts of ML trees for each genome
16S sequences, as a variability indicator.

rRNADif was created at Ivan Franko National University of Lviv.

Usage: -i|--input [FILE]; -d|--database_file [FILE]; -p|--project_name [NAME];
       [-n|--database_name [NAME]]; [-s|--sequence [FILE]]; [--plot][--tree]; 
       [--makedb_file [FILE]]; [--makedb_genomes [FOLDER]]; [--split_download [NUM]]
       [--step [1,2,3,4,0]]; [-m|--msa_alg [mafft, muscle,clustalo]]; 
       [-t|--tree_alg [fasttree, iqtree]]; [--only_phylo]; [-h|--help]; 
       [-v|--version]

The required arguments are: -i ; -d; -p.
If you are providing MSA or nwk file, please condider using -s flag to 
reconstruct phylogeny

    -i|--input        - input file in fna format. Genome sequence of interest.

    -d|--database_file- file to subset master database for branch length values 
                        extraction or path to custom database. A csv file, 
                        downloaded from "Browse by organism" component of NCBI 
                        Genome. For more information please see README at GitHub 
                        page. 

    -p|--project_name - prefix to name results files and folder

    -n|--database_name- name of database to use. Default: Actinobacteria. 
                        Databases can be found at datasets/computed folder. 
                        Folder name is database name.

    -s|--sequence     - consider providing sequence files for phylogeny run, if 
                        step greater than 0. Sequences must be in fasta format 
                        and have the same names as in provided .nwk or .fasta 
                        (MSA in fasta format) files.

    --plot            - save density plot of non-outlier values of used database

    --tree            - make phylogenetic tree for analyzed 16S rRNA sequences, 
                        where non-outlier organisms are represented with 1 
                        sequence (organism representative) and outliers - with 
                        all 16S sequences. Default algorithms fasttree and mafft.

    --redo            - force use provided directory (based on project name). 
                        All files will be rewritten and the pipeline will be 
                        rerun.

    --makedb_file     - make custom database with csv file from "Browse by 
                        genome". For more information look at GitHub page.
    
    --split_download  - with provided makedb_file, downloads genomes in split, 
                        preserving space. After 16S rRNA annotation is done, 
                        deletes downloaded genomes and begins a new cycle. Must 
                        go	with numerical value after, of how many genomes are 
                        preffered to be downloaded in one batch.

    --makedb_genomes  - make custom database in provided folder with genomes' 
                        sequences. Note: filenames will be treated as organism 
                        names!!!!

    --step: 1     - Run only barrnap on provided sequence file. 
            2     - Begin from provided multi-fasta file with 16S rRNA sequences.
            3     - MSA is provided via --input flag (fasta format).
            4     - Read nwk file and extract branch length for subsequent 
                    analysis
            0     - Run all (default)

Examples: sh rnadif.sh -i Genome.fna -d Streptomyces.csv -p Streptomyces --step 
          2 (this will omit barrnap run)

          sh rnadif.sh -i Genome.fna -d Streptomyces.csv -p Streptomyces --step 
          3 (this will omit MSA run)

          sh rnadif.sh -i Genome.fna -d Streptomyces.csv -p Streptomyces --step 
          1 (this will only run barrnap)		  

    -m|--msa_alg      - choose MSA algorithm to use in the pipeline. Can be also 
                        used for database computing
                        Options:
                        mafft (default)
                        muscle (need to be manually installed
                        clustalo (need to be manually installed)

    -t|--tree_alg     - choose phylogenetic inference program to run. Can be 
                        also used for database computing
                        Options:
                        fasttree (default)
                        iqtree (need to be manually installed. Works only with 
                        >=3 sequences in a file)
                        raxml (need to be manually installed)

    --only_phylo      - use chosen tree and msa tool only to construct final 
                        phylogenetic tree. For 16S tree fasttree and mafft will 
                        be used. 

    -h|--help - print this message and exit

    -v|--version - print version and exit. 


GitHub page: https://github.com/pavlohrab/rRNADif

Feel free to post any issues!


EOF
}

# Parsing script arguments
while (( "$#" )); do
  case "$1" in
    -i|--input)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            INPUT=$2
            INPUT_BOOL=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -d|--database_file)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            DATABASE=$2
            DATABASE_BOOL=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -p|--project_name)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            NAME=$2
            NAME_BOOL=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -n|--database_name)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            DATABASE_FOLDER=$2
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -s|--sequence)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            SEQUENCE=$2
            SEQUENCE_B=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    --plot)
      PLOT=1
      shift
      ;;
    --redo)
      REDO=1
      shift
      ;;
    --tree)
      TREE=1
      shift
      ;;
    --step)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            STEP=$2
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    --makedb_file)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            MAKEDB=$2
            MAKEDB_BOOL=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    --split_download)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            SPLIT=$2
            SPLIT_BOOL=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            SPLIT=1000
            echo "Using 1000 as an argument"
            shift
          fi
      ;;
    --makedb_genomes)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            MAKEDB_FOLDER=$2
            MAKEDB_BOOL_1=1
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -m|--msa_alg )
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            MSA_ALG=$2
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    -t|--tree_alg)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            TREE_ALG=$2
            shift 2
          else
            echo "Error: Argument for $1 is missing" >&2
            exit 1
          fi
      ;;
    --only_phylo)
      PHYLO=1
      shift
      ;;
    -h|--help)
      help
      exit 1
      shift
      ;;
    -v|--version)
      echo "Version number is: ${VERSION}"
      exit 1
      shift
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      help
      exit 1
      ;;
  esac
done

# Check if squence files are passed. If no, then check if the step == 2 (then 
# the input is sequences). If the step is != 0 and phylogeny computation is on
# then throw an error 
if  [ "$SEQUENCE_B" -eq 0 ]; then
    # de facto is step 0. Sequences can be used for tree computation
    if [ "$STEP" -eq 2 ]; then
        SEQUENCE_B=1
    elif
    # Throw an error. Need to provide sequences
        [ "$TREE" -eq 1 ] && [ "$STEP" -gt 0 ]; then
        echo ""
        echo "Flags --tree and step != 0 cannot be used in one run"
        echo "Please consider providing sequences with -s flag"
        echo "Sequence names must be identical!"
        hr
        exit 1
    fi
fi

# If no database is intented to be made (boolean zeros) then make sure that 
# alll the mandatory arguments are passed to the script
if [ "$MAKEDB_BOOL" -eq 0 ] && [ "$MAKEDB_BOOL_1" -eq 0 ]; then
    if [[ "$NAME_BOOL" -eq 0 ]] ; then
        echo "Please provide -p argument (project name)."
        echo "For more information use -h flag"
        hr
        exit 1
    fi
    if [[ "$INPUT_BOOL" -eq 0 ]] ; then
        echo "Missing input genome sequence file. Please use -i flag"
        echo "For more information use -h flag"
        hr
        exit 1
    fi
    if [[ "$DATABASE_BOOL" -eq 0 ]] ; then
        echo "No database is provided. Please use -d or -cd flag."
        echo "For more information use -h flag"
        hr
        exit 1
    fi
fi

# If user is generating a database, then print the below message
if [[ "$MAKEDB_BOOL" -eq 1 ]] || [[ "$MAKEDB_BOOL_1" -eq 1 ]]; then
    echo "Only making custom database. All other arguments are ignored"
    echo "To run the program with custom database please use -d flag next time"
fi

# Go into database creation. 
if [ "$MAKEDB_BOOL" -eq 1 ] || [ "$MAKEDB_BOOL_1" -eq 1 ]; then
    # Check if we are making database from file. If so, go into this code
    if [ "$MAKEDB_BOOL" -eq 1 ] ; then 
        echo 'Downloading and making database'
        # Make tmp folder and copy input there
        mkdir tmp_database 
        cp $MAKEDB tmp_database
        # Navigate to that directory
        cd tmp_database
        # The file itself is the last bit of path string. Extract that.
        MAKEDB_NAME=$(echo "$MAKEDB" | awk -F '/' '{print $NF}')		
        hr
        # Some logic about if the phylogeny and msa algorithms are passed through
        # arguments and if split download is needed. 

        # TO-DO: Do we need this logic? The MSA_ALG and TREE_ALG variables are 
        # set to defaults.....
        if [ "$PHYLO" -eq 0 ]; then
            if [ "$SPLIT_BOOL" -eq 1 ]; then
                sh ../scripts/annotate_and_rename_sequences.sh $MAKEDB_NAME \
                ../scripts/ $MSA_ALG $TREE_ALG $SPLIT
            else
                sh ../scripts/annotate_and_rename_sequences.sh $MAKEDB_NAME \
                ../scripts/ $MSA_ALG $TREE_ALG 0
            fi
        else
            if [ "$SPLIT_BOOL" -eq 1 ]; then
                sh ../scripts/annotate_and_rename_sequences.sh $MAKEDB_NAME \
                ../scripts/ mafft fasttree $SPLIT
            else
                sh ../scripts/annotate_and_rename_sequences.sh $MAKEDB_NAME \
                ../scripts/ mafft fasttree 0
            fi
        fi
        # Get the actual filename 
        FILENAME="$(basename -s .csv $MAKEDB_NAME)"
        # Make directory in datasets/computed with this filename and move-rename
        # all the generated result files there.
        mkdir "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_all.csv"  "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_one.csv"  "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_mean_outliers.csv"  "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_median_outliers.csv"  "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_no_ouliers.csv"  "../datasets/computed/${FILENAME}"
        mv "${FILENAME}_16S_rna" "../datasets/computed/${FILENAME}/rrna"
        mv "${FILENAME}_one_rna" "../datasets/computed/${FILENAME}/one_rrna"
        mv "${FILENAME}_no_rrna" "../datasets/computed/${FILENAME}/no_rrna"
        # Remove all tmp folder
        cd .. && rm -rf mkdir tmp_database 
        # Remain csv file which was used to make database
        cp $MAKEDB datasets/downloaded_csv/$MAKEDB
        hr
        echo "Database is created and accesible via -n flag"
        echo "For database name please see datasets/computed/folder-name"
        hr
        exit 1
    elif [ "$MAKEDB_BOOL_1" -eq 1 ] ; then 
        echo "Making database from sequence folder"
        # Make tmp folder and copy inputs there
        mkdir tmp_database
        cp "${MAKEDB_FOLDER}"/* tmp_database
        MAKEDB_NAME=$(basename "$MAKEDB_FOLDER")
        # Make the folder in datasets/computed
        mkdir "datasets/computed/${MAKEDB_NAME}"
        cd tmp_database
        ls * > list.txt
        cp ../scripts/annotate_and_rename_sequences_folder.sh ./
        # TO-DO: Do we need this?
        if [ "$PHYLO" -eq 0 ]; then
            sh annotate_and_rename_sequences_folder.sh $MAKEDB_NAME ../scripts/ \
            $MSA_ALG $TREE_ALG
        else
            sh annotate_and_rename_sequences_folder.sh $MAKEDB_NAME ../scripts/ \
            mafft fasttree
        fi
        # Move-rename the result files into datasets/computed folder
        mv "${MAKEDB_NAME}_all.csv"  "../datasets/computed/${MAKEDB_NAME}"
        mv "${MAKEDB_NAME}_one.csv"  "../datasets/computed/${MAKEDB_NAME}"
        mv "${MAKEDB_NAME}_mean_outliers.csv"  "../datasets/computed/${MAKEDB_NAME}"
        mv "${MAKEDB_NAME}_median_outliers.csv"  "../datasets/computed/${MAKEDB_NAME}"
        mv "${MAKEDB_NAME}_no_ouliers.csv"  "../datasets/computed/${MAKEDB_NAME}"
        mv "${MAKEDB_NAME}_16S_rna" "../datasets/computed/${MAKEDB_NAME}/rrna"
        mv "${MAKEDB_NAME}_one_rrna" "../datasets/computed/${MAKEDB_NAME}/one_rrna"
        mv "${MAKEDB_NAME}_no_rrna" "../datasets/computed/${MAKEDB_NAME}/no_rrna"
        # rm tmp folder
        cd .. && rm -rf mkdir tmp_database 
        # Copy csv file. used for database generation into datasets/downloaded_csv 
        cp $MAKEDB datasets/downloaded_csv/$MAKEDB
        hr
        echo "Database is created and accesible via -n flag"
        echo "For database name please see datasets/computed/folder-name"
        hr
        exit 1
    fi
fi

# The main script if no database creation is need to be done

# Get the database and input names. First get the files is 
DATABASE_NAME=$(echo "$DATABASE" | awk -F '/' '{print $NF}')
DATABASE_NAME_2=$(basename -s .csv $DATABASE_NAME)
INPUT_NAME=$(echo "$INPUT" | awk -F '/' '{print $NF}')
INPUT_NAME_2=$(basename "$INPUT_NAME" | sed 's/\(.*\)\..*/\1/')

# Go into here if no database creation arguments were passed. No actual need
# for this check due to sequential nature of code execution.
if [ "$MAKEDB_BOOL" -eq 0 ] && [ "$MAKEDB_BOOL_1" -eq 0 ]; then
    # Firstly check if folder withe same project name exists. If so stop the 
    # execution and print error message
    if [ -d "${NAME}_results" ]; then
        if [ "$REDO" -eq 0 ]; then
            echo "Results directory with provided name already exists"
            echo "Please change project name or consider using --redo flag"
            hr
            exit 1
        # If the user passed --redo argument, rewrite result files in that folder
        elif [ "$REDO" -eq 1 ]; then
            echo "Rewriting files in ${NAME}_results"
            rm -rf "${NAME}_results"
            mkdir "${NAME}_results" 
            mkdir "${NAME}_results/database_rrnas"
            mkdir "${NAME}_results/input_files"
            mkdir "${NAME}_results/results"
        fi
    fi 

    # If no folder with project name was found, create one.
    if [ ! -d "${NAME}_results" ]; then
        mkdir "${NAME}_results" 
        mkdir "${NAME}_results/database_rrnas"
        mkdir "${NAME}_results/input_files"
        mkdir "${NAME}_results/results"
    fi

    # Copy input flies into project folder 
    cp "$INPUT" "${NAME}_results/input_files" && cp "$DATABASE" "${NAME}_results/input_files"

    # Preprocess the input database file. Rename the sequences and search 
    # for them in the precomputed database.
    python scripts/dataframe_preprocess_for_subset.py \
    "${NAME}_results/input_files/${DATABASE_NAME}" \
    "${NAME}_results/results/${DATABASE_NAME_2}_used_organisms.csv" \
    "datasets/computed/${DATABASE_FOLDER}" $DATABASE_FOLDER \
    "${NAME}_results/results/${DATABASE_NAME_2}_one.csv"
    
    # If returned dataframe is empty then no sequences are in the precomputed 
    # database
    if [ "$(wc -l  "${NAME}_results/results/${DATABASE_NAME_2}_used_organisms.csv" | awk '{print $1}')" -eq 1 ]; then
        echo "No organisms found in chosen database with ${DATABASE_NAME} input"
        echo "Please condider choosing different database with -n flag"
        hr
        exit 1
    fi

    # Recompute statistics for generated subset with help of R script
    echo "Creating database subset with provided csv file"
    Rscript --vanilla scripts/statistics_for_subset.R \
    $DATABASE_NAME_2 \
    "${NAME}_results/results/${DATABASE_NAME_2}_used_organisms.csv" "$(pwd)"  \
    "${NAME}_results/results"
    hr

    # Do some preprocessing if phylogenetic tree should be built. 
    # Actually just copy the fasta files using the subset from master dataframe 
    if [ "$TREE" -eq 1 ]; then
        echo ""

        echo "Preparing rrnas for phylogenetic inference..."
        echo "This can take some time"

        # Actually generates txt files where organism names are listed in every 
        # line. For non outlier, outliers and single sequence different txt files 
        # are generated
        python scripts/preprocess_dataframes.py \
        "${NAME}_results/results/${DATABASE_NAME_2}" \
        "${NAME}_results/results" "$(pwd)"
        
        # Make a folder for filtered 16S rRNAs 
        mkdir "${NAME}_results/database_rrnas/cleaned"

        # Copy fasta files for those organisms who need to be filtered (non 
        # outliers) into dayabase_rrnas folder. 
        while read p; do cp "datasets/computed/${DATABASE_FOLDER}/rrna/$p.fasta" \
        "${NAME}_results/database_rrnas"; done < "${NAME}_results/results/remain_one.txt" 
        
        # Copy outlier rrnas and single files with rrnas into cleaned folder.
        # We are goind to use them all,
        while read p; do cp "datasets/computed/${DATABASE_FOLDER}/rrna/$p.fasta" \
        "${NAME}_results/database_rrnas/cleaned"; done < "${NAME}_results/results/outliers_mean.txt"
        while read p; do cp "datasets/computed/${DATABASE_FOLDER}/rrna/$p.fasta" \
        "${NAME}_results/database_rrnas/cleaned"; done < "${NAME}_results/results/outliers_median.txt"
        while read p; do cp "datasets/computed/${DATABASE_FOLDER}/one_rrna/$p.fasta" \
        "${NAME}_results/database_rrnas/cleaned"; done < "${NAME}_results/results/only_one.txt"

        # Add number to a fasta header, if it it dukicate in this file. Because 
        # all sequence in a file are named by organism name, then this step is 
        # necessary to run phylogeny
        cd "${NAME}_results/database_rrnas/cleaned" && for i in *.fasta; do perl -pe \
        's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' $i > $i.clean; done && cd ../../../

        # Rewrite fastas for non outlier organisms to remain only one sequence
        cat "${NAME}_results/results/remain_one.txt" | parallel python scripts/remain_one.py {} \
        "${NAME}_results/database_rrnas/cleaned" "${NAME}_results/database_rrnas"

        # Remove innitial fasta files
        rm -f "${NAME}_results"/database_rrnas/*.fasta
        rm -f "${NAME}_results"/database_rrnas/cleaned/*.fasta
        hr
    fi

    # Check the step variable (the --step argument passed) to know from where
    # to run an analysis. If step == 0 (default), then all the below chunks of 
    # code would run (all the analysis pipeline) 
    # If step == 0 or 1, then run barrnap
    if [ "$STEP" -eq 0 ] || [ "$STEP" -eq 1 ]; then
        echo ""
        echo "Step 1: Annotating genome sequence with barrnap..."
        # Barnap run
        barrnap --quiet  --outseq "${NAME}_results/results/${INPUT_NAME_2}_all_rrna.fasta" \
        "${NAME}_results/input_files/${INPUT_NAME}"

        # Linearize sequences, extract only 16S 
        cat "${NAME}_results/results/${INPUT_NAME_2}_all_rrna.fasta" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  \
        | grep -A 1 "^>16S" | sed 's/:/_/g' > "${NAME}_results/results/${INPUT_NAME_2}_16S_rrna.fasta"

        # Rename sequences and remove duplicate names if more than one 16S is
        # present (add number there)
        cd "${NAME}_results/results" && awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-6); next} 1' "${INPUT_NAME_2}_16S_rrna.fasta" \
        | sed 's/_16S_rrna//g' | perl -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; '\
         > "${INPUT_NAME_2}_16S_clean_rrna.fasta" && cd ../../
        # If step was passed as 1, then the analysis is complete. Exit the program
        if [ "$STEP" -eq 1 ]; then
            echo ""
            echo "16S rRNA sequences are annotated!"
            echo "For full analysis please run with --step 0"
            hr
            exit 1
        fi
    fi
    hr

    # Check the step. If 0 or 2 then run
    if [ "$STEP" -eq 0 ] || [ "$STEP" -eq 2 ]; then
        # Check if step ==2. If so, then some preparation is needed 
        if  [ "$STEP" -eq 2 ] ; then 
            echo ""
            echo "Step 2: Reading input 16S rRNA files..."

            # Check the algoritms which the user prefers and run MSA
            # The output files are the same to make their furhter call easier
            if [ "$MSA_ALG" = clustalo ]; then
                clustalo -i "${NAME}_results/input_files/${INPUT}" \
                 -o "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            elif [ "$MSA_ALG" = muscle ]; then
                muscle -quiet  = "${NAME}_results/input_files/${INPUT}" \
                 -out "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            elif [ "$MSA_ALG" = mafft ]; then
                mafft --quiet "${NAME}_results/input_files/${INPUT}" > \
                "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            else
                echo "Named MSA algorithm was not found! "
                echo "Please check the spelling and look in help for available options"
                hr
                exit 1
            fi
            # Copy the input sequences to use in the phylogeny 
            if [ "$SEQUENCE_B" -eq 1 ]; then
                cp  "${NAME}_results/input_files/${INPUT}" \
                "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta"
            fi
            # Return STEP variable to 0 (to make subsequent analysis run)
            STEP=0
        else
            echo ""
            echo "Step 2: Running ${MSA_ALG} for MSA generation"
            
            # Running MSA for annotated 16S sequences with chosen program  
            if [ "$MSA_ALG" = clustalo ]; then
                clustalo -i "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta"\
                  -o "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            elif [ "$MSA_ALG" = muscle ]; then
                muscle -quiet -in "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta" \
                 -out "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            elif [ "$MSA_ALG" = mafft ]; then
                mafft --quiet "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta"\
                 > "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            else
                echo "Named MSA algorithm was not found! "
                echo "Please check the spelling and look in help for available options"
                hr
                exit 1
            fi
        fi
    fi
    hr

    # Check the step is 0 or 3 then run the code below
    if  [ "$STEP" -eq 0 ] || [ "$STEP" -eq 3 ]; then
        # If the step ==3, thsn read the input and run phylogeny
        if [ "$STEP" -eq 3 ]; then 
            echo ""
            echo "Step 3: Using provided MSA file..."
            # Run chosen phylogeny program with input MSA
            # Also rename the output, if exention is not .nwk
            if  [ "$TREE_ALG" = fasttree ]; then
                fasttree -quiet -nt -gtr -out \
                "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk" \
                "${NAME}_results/input_files/${INPUT}"
            elif [ "$TREE_ALG" = iqtree ]; then
                iqtree --quiet -T 1 -m GTR -s "${NAME}_results/input_files/${INPUT}"
                cp "${INPUT}"*.treefile "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk" 
            elif [ "$TREE_ALG" = raxml ]; then
                raxmlHPC -s "${NAME}_results/input_files/${INPUT}" -n \
                "${INPUT_NAME_2}".tmp -m  GTRCAT --print-identical-sequences
                cp RAxML_bestTree.*"${INPUT_NAME_2}".tmp \
                "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk"
                rm *.tmp
            else
                echo "Named TREE algorithm was not found! "
                echo "Please check the spelling and look in help for available options"
                hr
                exit 1
            fi
            # Copy provided 16S sequences for latter use
            if [ "$SEQUENCE_B" -eq 1 ]; then
                cp  "${SEQUENCE}" "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta"
            fi
            # Make STEP variable ==0 to run subsequent analysis
            STEP=0
        else
            echo ""
            echo "Step 3: Computing phylogenetic tree with ${TREE_ALG}..."
            # Run chosen phylogeny program for generated MSA
            # Also rename the output, if exention is not .nwk
            if  [ "$TREE_ALG" = fasttree ]; then
                fasttree -quiet -nt -gtr -out "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk" \
                "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
            elif [ "$TREE_ALG" = iqtree ]; then
                iqtree --quiet -T 1 -m GTR -s "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta"
                cp "${INPUT_NAME_2}"*.treefile "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk"
            elif [ "$TREE_ALG" = raxml ]; then
                raxmlHPC -s "${NAME}_results/results/${INPUT_NAME_2}_16S_mafft.fasta" -n \
                "${INPUT_NAME_2}".tmp -m  GTRCAT --print-identical-sequences
                cp RAxML_bestTree.*"${INPUT_NAME_2}".tmp "${NAME}_results/results/${INPUT_NAME_2}_16S.nwk"
                rm *.tmp
            else
                echo "Named TREE algorithm was not found! "
                echo "Please check the spelling and look in help for available options"
                hr
                exit 1
            fi
        fi
    fi
    hr

    # If step ==0 or 4 run the below code
    if  [ "$STEP" -eq 0 ] || [ "$STEP" -eq 4 ]; then
        if [ "$STEP" -eq 4 ]; then 
            echo ""
            echo "Step 4: Extracting branch length from provided nwk file..."
            # Rename the input file to have .nwk exention
            cp "${NAME}_results/input_files/${INPUT}"  "${NAME}_results/results/${INPUT_NAME_2}.nwk"
            # Extract branch length for and compute mean/median for the provided tree 
            Rscript --vanilla scripts/compute_dataframes_1_seq.R $INPUT_NAME_2 "${NAME}_results/results"
            if [ "$SEQUENCE_B" -eq 1 ]; then
                cp  "${SEQUENCE}" "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta"
            fi
            STEP=0
        else
            echo ""
            echo "Step 4: Extracting branch length from file."
            Rscript --vanilla scripts/compute_dataframes_1_seq.R $INPUT_NAME_2 "${NAME}_results/results"
        fi
    fi
    hr
    echo "Combining data from database and provided sequence"
    Rscript --vanilla scripts/combine_dataframes.R "${INPUT_NAME_2}_all.csv" \
    "${DATABASE_NAME_2}_all.csv" "${NAME}_results/results"

    if [ "$TREE" -eq 1 ]; then
        echo ""
        echo "Running phylogenetic inference..."
        mkdir "${NAME}_results/results/phylogeny"
        cp "${NAME}_results/results/${INPUT_NAME_2}_16S_clean_rrna.fasta" \
        "${NAME}_results/database_rrnas/cleaned/${INPUT_NAME_2}_16S_clean_rrna.clean"
        cat "${NAME}_results"/database_rrnas/cleaned/*.clean > \
        "${NAME}_results/results/phylogeny/FINAL.fasta"
        if [ "$MSA_ALG" = clustalo ]
            then
            clustalo -i "${NAME}_results/results/phylogeny/FINAL.fasta" \
             -o "${NAME}_results/results/phylogeny/FINAL.mafft"
        elif [ "$MSA_ALG" = muscle ]
            then
            muscle  -in "${NAME}_results/results/phylogeny/FINAL.fasta"\
             -out "${NAME}_results/results/phylogeny/FINAL.mafft"
        elif [ "$MSA_ALG" = mafft ]
            then
            mafft --quiet "${NAME}_results/results/phylogeny/FINAL.fasta"\
            > "${NAME}_results/results/phylogeny/FINAL.mafft"
        else
            echo "Named MSA algorithm was not found! "
            echo "Please check the spelling and look in help for available options"
            hr
            exit 1
        fi

        if  [ "$TREE_ALG" = fasttree ]
            then
            fasttree -quiet -nt -gtr -out "${NAME}_results/results/phylogeny/FINAL.nwk"\
             "${NAME}_results/results/phylogeny/FINAL.mafft"
        elif [ "$TREE_ALG" = iqtree ]
            then
            iqtree  -T 1 -m GTR -s "${NAME}_results/results/phylogeny/FINAL.mafft"
            cp FINAL.mafft.treefile "${NAME}_results/results/phylogeny/FINAL.nwk"
        elif [ "$TREE_ALG" = raxml ]
            then
            raxmlHPC -s "${NAME}_results/results/phylogeny/FINAL.mafft" -n \
            FINAL.tmp -m  GTRCAT --print-identical-sequences
            cp RAxML_bestTree.FINAL.tmp "${NAME}_results/results/phylogeny/FINAL.nwk"
        else
            echo "Named TREE algorithm was not found! "
            echo "Please check the spelling and look in help for available options"
            hr
            exit 1
        fi
    fi

    if [ "$PLOT" -eq 1 ]; then
        echo ""
        echo "Generating Rplot for non_outlier values..."
        Rscript --vanilla scripts/plot.R Results_no_ouliers.csv "${NAME}_results/results"
    fi
    hr
    cd "${NAME}_results/results"
    mkdir "database_used_data"
    mkdir "input_sequence_individual results"
    mv $DATABASE_NAME_2* "database_used_data"
    mv $INPUT_NAME_2* "input_sequence_individual results"
    if [ "$TREE" -eq 1 ]; then
        rm  *.txt
    fi
    rm 	Rplots.pdf
    RESULT_NAME="${INPUT_NAME_2}"
    RESULT_MEAN=$(cat Results_all.csv | grep "${INPUT_NAME_2}" | awk -F ',' '{print $2}')
    RESULT_MEDIAN=$(cat Results_all.csv | grep "${INPUT_NAME_2}" | awk -F ',' '{print $3}')
    echo ""
    echo "Results are:"
    echo "Species name: ${RESULT_NAME}"
    echo "Mean value: ${RESULT_MEAN}"
    echo "Meadian value: ${RESULT_MEDIAN}"
    if [ "$(cat Results_no_ouliers.csv | grep "${INPUT_NAME_2}")" ]; then
        echo "The results are within no_outlier values"
    fi
    if [ "$(cat Results_mean_outliers.csv | grep "${INPUT_NAME_2}")" ]; then
        echo "The results are within mean outlier values"
    fi
    if [ "$(cat Results_median_outliers.csv | grep "${INPUT_NAME_2}")" ]; then
        echo "The results are within meadian values"
    fi
    cd ../../
    rm 	Rplots.pdf
    hr
    echo ""
    echo "The analysis is complete. Results are in the project-name/results folder."
    echo ""
fi