#!/bin/bash

#command
usage() {
    echo "Usage: $0 --working-dir <WORKING_DIR> --input <INPUT_FILES_DIR> --accession <ACCESSION> --kmer-size <KMER_SIZE> --organism <ORGANISM> [--nowget --reads-dir <READS-DIR>]"
    exit 1
}

#default empty values for wget flag
OPTIONAL_FLAG=""
READS_DIR=""

#command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --working-dir)
            WORKING_DIR=$2
            shift 2
            ;;
        --input)
            INPUT=$2
            shift 2
            ;;
        --accession)
            ACCESSION=$2
            shift 2
            ;;
        --kmer-size)
            KMER_SIZE=$2
            shift 2
            ;;
        --organism)
            ORGANISM=$2
            shift 2
            ;;
        --no-wget)
            OPTIONAL_FLAG="--no-wget"
            shift
            ;;
        --reads-dir)
            READS_DIR=$2
            shift 2
            ;;
        *)
            echo "unknown flag: $1"
            usage
            ;;
    esac
done

#check for required arguments
if [[ -z "$WORKING_DIR" || -z "$INPUT" || -z "$ACCESSION" || -z "$KMER_SIZE" || -z "$ORGANISM" ]]; then
    usage
fi

#make sure --reads-dir is specified if --no-wget is included
if [[ "$OPTIONAL_FLAG" == "--no-wget" && -z "$READS_DIR" ]]; then
    echo "Error: --no-wget requires specifying --reads-dir."
    exit 1
fi

#make pipeline directory
DIR="$WORKING_DIR/pipelineproject_claudia_krueger"
mkdir -p "$DIR"

#find the downloaded scripts
USER_HOME=$(eval echo ~$USER)
SCRIPT_DIR=$(find "$USER_HOME" -type d -name "biopython_pipeline" -print -quit)

#determine if wget should be used to retrieve fasta files
if [[ "$OPTIONAL_FLAG" == "--no-wget" ]]; then
    mkdir -p "$DIR/reads"
else
    #if --no-wget not included, run the first Python script to retrieve files
    mkdir -p "$DIR/reads"
    python "$SCRIPT_DIR/scripts/01file_retrieving.py" --input "$INPUT" --dir "$DIR"
fi

#run the pipeline with the specified directory and accession

#find cds
python "$SCRIPT_DIR/scripts/02find_cds.py" --dir "$DIR" --accession "$ACCESSION"

#kallisto
if [[ "$OPTIONAL_FLAG" == "--no-wget" ]]; then
    python "$SCRIPT_DIR/scripts/03kallisto.py" --input "$INPUT" --dir "$DIR" --path "$READS_DIR" --no-wget
else
    python "$SCRIPT_DIR/scripts/03kallisto.py" --input "$INPUT" --dir "$DIR"
fi

#sleuth
python "$SCRIPT_DIR/scripts/04sleuth_input.py" --input "$INPUT" --dir "$DIR"

#sleuth
Rscript "$SCRIPT_DIR/scripts/04sleuth.R" --dir "$DIR"

#bowtie
python "$SCRIPT_DIR/scripts/05bowtie2.py" --input "$INPUT" --dir "$DIR" --accession "$ACCESSION" --path "$READS_DIR" --no-wget

#spades
python "$SCRIPT_DIR/scripts/06spades.py" --input "$INPUT" --dir "$DIR" --kmer "$KMER_SIZE"

#blast
python "$SCRIPT_DIR/scripts/07blast.py" --input "$INPUT" --dir "$DIR" --organism "$ORGANISM"