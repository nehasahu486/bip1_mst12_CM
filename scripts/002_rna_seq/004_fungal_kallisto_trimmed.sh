#!/bin/bash
#SBATCH -p tsl-medium                           # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH --wckey=talbot_core                     # project code
#SBATCH --mem 200G                              # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,fail              # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address



source package /tsl/software/testing/bin/kallisto-0.46.2


## STEP 1 - Index the transcriptome with Kallisto
kallisto index -i data/processed/MG8_cdna.idx /data_raw/Magnaporthe_oryzae.MG8.cdna.all.fa

## STEP 2 - Perform pseudoalignment using Kallisto on trimmed FASTQ files
mkdir -p kallisto_output_fungal

for file in trimmed/*_paired_1.fastq.gz; do
kallisto quant -i data/processed/MG8_cdna.idx \
 -o kallisto_output_fungal/$(basename "$file" _paired_1.fastq.gz) \
 --threads 24 \
 -b 100 \
 "$file" "${file/_paired_1.gz/_paired_2.gz}"
done