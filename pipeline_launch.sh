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

function align_sort_index() {
    # align with bwa mem, sort and index with samtools
    REF=$1
    INP=$2
    OUT=${INP%.*}.bam
    set -o pipefail
    return $(bwa mem $REF $INP | samtools sort > $OUT && samtools index $OUT)
}

# Arguments
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <brass_input.bedpe> <reference.fasta> <bam_file1> [<bam_file2> ...]"
    exit
fi

THIS_SCRIPT=$(abspath $0)
THIS_DIR=$(dirname $THIS_SCRIPT)
WORKING_DIR=$(/bin/pwd)
CHECKPOINT=$WORKING_DIR/logs/CHECKPOINT
mkdir -p $WORKING_DIR/logs
touch $CHECKPOINT
BRASS_INPUT=$1
REFERENCE=$2
shift
shift
BAMS=$@
echo
echo "############################"
echo "Breakpoint Assembly Pipeline"
echo $THIS_SCRIPT
echo
echo "############################"
echo
echo


OUTPUT=sv_predictions.txt
PIPELINE=$THIS_DIR
INPUT="${WORKING_DIR}/tigra_input.txt"

## STEP 1
# Convert BRASS input to TIGRA input
if [ ! -s $INPUT ] || ! grep -q "1_create_input" $CHECKPOINT; then
    echo "Writing TIGRA input file to $(/bin/pwd) $INPUT"
    echo
    if python ${PIPELINE}/convert_input.py ${BRASS_INPUT} ${INPUT}; then
        echo "Success"
        echo "1_create_input" >> $CHECKPOINT
        echo
    else
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Failed to create input file"
    fi
fi

## STEP 2
# Run TIGRA
if ! grep -q "2_tigra" $CHECKPOINT; then

    echo "Running TIGRA. Contigs will be written to ${WORKING_DIR}/contigs.fa"
    echo "Reads will be written to ${WORKING_DIR}/reads_fasta/"
    echo "Tigra logs will be written to ${WORKING_DIR}/logs/tigra.log"
    mkdir -p ${WORKING_DIR}/reads_fasta
    mkdir -p ${WORKING_DIR}/logs
    # NB: TIGRA step will take 2-3 times longer with the -m flag set because it has to search
    # for mate paired reads. In tests this means ~3 minutes per SV with -m set, and ~1 minute otherwise.
    if ! command -v tigra-sv > /dev/null; then
        echo "tigra-sv is not in the PATH, pipeline aborted"
        exit 2
    fi
    if tigra-sv -l 1000 \
                -a 100 \
                -o ${WORKING_DIR}/contigs.fa \
                -R ${REFERENCE} \
                -w 400 \
                -p 40000 \
                -h 1000 \
                -d -I reads_fasta \
                -b -m -g 2 ${INPUT} ${BAMS} 2> ${WORKING_DIR}/logs/tigra.log; then
        echo "Success"
        echo "2_tigra" >> ${CHECKPOINT}
        echo
    else
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Failed to assemble contigs with TIGRA"
        exit 2;
    fi
fi

## STEP 2a
# Check that contig names are valid
if ! $(python -c "import skbio" 2&1> /dev/null); then
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "Could not load library skbio"
    exit 1
fi
echo
echo "Validating contig names"
python $PIPELINE/validate_contig_names.py contigs.fa tmp full_contig_names.txt && mv tmp contigs.fa
echo "Success"
echo

## STEP 3
# Align new contigs back to the reference
if ! grep -q "3_align_contigs" $CHECKPOINT; then
    echo "Aligning contigs to reference with BWA-MEM"
    if align_sort_index ${REFERENCE} ${WORKING_DIR}/contigs.fa; then
        echo "Success"
        echo "3_align_contigs" >> ${CHECKPOINT}
        echo
    else
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Failed to align contigs"
        exit 3
    fi
fi

## STEP 4
# Align reads to the reference
if ! grep -q "4_align_reads" $CHECKPOINT; then
    echo "Aligning reads to reference with BWA-MEM"
    cat ${WORKING_DIR}/reads_fasta/*.fa > ${WORKING_DIR}/reads_fasta/allreads.fa
    if align_sort_index ${REFERENCE} ${WORKING_DIR}/reads_fasta/allreads.fa; then
        echo "Success"
        echo "4_align_reads" >> $CHECKPOINT
        echo
    else
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!"
        echo "Failed to align reads"
        exit 4
    fi
fi

## STEP 5
# Run python scripts
if ! grep -q "5_map_contigs" $CHECKPOINT; then
    if ! $(python -c "import pysam" 2&1> /dev/null); then
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "Could not load library pysam"
        exit 1
    fi

    python $PIPELINE/extract_supplementary_reads.py reads_fasta/allreads.bam reads_fasta/suppreads.bam
    python $PIPELINE/map_contigs.py contigs.bam contigs.fa > ${OUTPUT}
    if [ -s "$OUTPUT" ]; then
        echo "5_map_contigs" >> $CHECKPOINT
    else
        echo "Failed to map contigs"
        exit 5;
    fi
fi
echo "Done."
