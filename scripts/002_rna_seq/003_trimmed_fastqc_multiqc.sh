#!/bin/bash
#SBATCH -p tsl-medium                           # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH --wckey=talbot_core                     # project code
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,failL             # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address


## STEP 1 - Define input and output directories
input_dir="trimmed"
fastqc_output="trimmed_fastqc"
multiqc_output="trimmed_multiqc"

## STEP 2 - Run FastQC on each trimmed FASTQ file
source package /tsl/software/testing/bin/fastqc-0.11.5

mkdir -p "$fastqc_output"
for file in "$input_dir"/*.fastq.gz; do
    fastqc "$file" -o "$fastqc_output"
done

## STEP 3 - Run MultiQC on the FastQC output directory
source package /tgac/software/testing/bin/multiqc-1.5_ei
mkdir -p "$multiqc_output"
multiqc -s "$fastqc_output" -o "$multiqc_output"