#!/bin/bash

function abspath() {
    # generate absolute path from relative path
    # $1     : relative filename
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    fi
}

# Arguments
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <brass_input.bedpe> <reference.fasta> <bam_file1> [<bam_file2> ...]"
    exit
fi
BRASS_INPUT=$1
REFERENCE=$2
shift
shift
BAMS=$@
OUTPUT=sv_predictions.txt
PIPELINE="/lustre/scratch112/sanger/casm/dog_n_devil/kg8/projects/structvar/tigra/pipeline"
INPUT="tigra_input.txt"

# Convert BRASS input to TIGRA input
echo "Writing TIGRA input file to $(/bin/pwd) $INPUT"
python ${PIPELINE}/convert_input.py ${BRASS_INPUT} ${INPUT}

# Run TIGRA
echo "Running TIGRA. Contigs will be written to $(/bin/pwd)/contigs.fa"
echo "Reads will be written to $(/bin/pwd)/reads_fasta/"
echo "Tigra logs will be written to $(/bin/pwd)/logs/tigra.log"
mkdir -p reads_fasta
mkdir -p logs
# NB: TIGRA step will take 2-3 times longer with the -m flag set because it has to search
# for mate paired reads. In tests this means ~3 minutes per SV with -m set, and ~1 minute otherwise.
tigra-sv -l 1000 \
         -a 100 \
         -o contigs.fa \
         -R ${REFERENCE} \
         -w 400 \
         -p 40000 \
         -h 1000 \
         -d -I reads_fasta \
         -b -m -g 2 ${INPUT} ${BAMS} 2> logs/tigra.log

# Align new contigs back to the reference
echo "Aligning contigs to reference with BWA-MEM"
bwa mem ${REFERENCE} contigs.fa | samtools sort > contigs.bam
samtools index contigs.bam
cat reads_fasta/*.fa > reads_fasta/allreads.fa
bwa mem ${REFERENCE} reads_fasta/allreads.fa | samtools sort > reads_fasta/allreads.bam
samtools index reads_fasta/allreads.bam
python $PIPELINE/extract_supplementary_reads.py reads_fasta/allreads.bam reads_fasta/suppreads.bam
python $PIPELINE/map_contigs.py contigs.bam contigs.fa > ${OUTPUT}

echo "Done."
