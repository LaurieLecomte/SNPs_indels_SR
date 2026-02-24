#!/bin/bash

# First split sample list to align into different files with:
# cd 04_bam
# ls -1 *.fq.gz > ../all_samples_for_alignment.txt
# cd ..
# mkdir samples_split
# split -a 4 -l 1 -d all_samples_for_alignment.txt samples_split/samples_split.

## With GNU Parallel
# ls -1 samples_split/* | parallel -k -j 8 ./01_scripts/utils/bwa_mem_single-end.sh 8 {}

## With GNU Parallel on slurm
# ls -1 samples_split/* | parallel -k -j 10 srun -c 8 --mem 30G -p large --time 21-00:00 -J bwaMem -o 10-log_files/bwaMEMsplit_%j.log ls -1 samples_split/* | parallel -k -j 10 ./01_scripts/utils/bwa_mem_single-end.sh 8 {}

# Global variables
GENOMEFOLDER="03_genome"
GENOME="genome.fasta"
DATAFOLDER="04_bam"
NCPU="$1"
SAMPLE_FILE="$2"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

# Modules
#module load bwa
#module load samtools

# Index genome if not alread done
# bwa index "$GENOMEFOLDER"/"${GENOME%.fasta}"

cat "$SAMPLE_FILE" |
while read file
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename "$file")
    sample=$(basename "$name" .fastq.gz)
    ID="@RG\tID:"$sample"\tSM:"$sample"\tPL:Illumina"

    # Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$DATAFOLDER"/"${name%.fq.gz}".bam

    # Samtools sort
    samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam

    # Cleanup
    rm "$DATAFOLDER"/"${name%.fq.gz}".bam
done
