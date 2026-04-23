#!/bin/bash
#SBATCH -p tsl-medium                           # queue
#SBATCH --wckey=talbot_core                     # project_code
#SBATCH -N 1                                    # number of nodes
#SBATCH --cpus-per-task=24                      # number of cores per task
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,failL             # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address



## Run FastQC
source package /tsl/software/testing/bin/fastqc-0.11.5

# Define input directory
input_dir1="/tsl/data/reads/ntalbot/timecourse_rnaseq_bip1_mst12/"
input_dir2="/tsl/data/reads/ntalbot/guy11_rnaseq_appressorium" #From Neftaly's Vts1 experiments

# Define output directory for fastqc
fastqc_output="fastqc_output"

# Create output directory
mkdir -p "$fastqc_output"

# Run FastQC on each FASTQ file found in both input directories and their subdirectories
find "$input_dir1" -name "*.fq.gz" | while read file; do
    fastqc "$file" -o "$fastqc_output"
done

## Run MultiQC
source package /tgac/software/testing/bin/multiqc-1.5_ei

# Define output directories for multiqc
multiqc_output="multiqc_output"

# Create output directory
mkdir -p "$multiqc_output"

# Run MultiQC on the FastQC output directory
multiqc -s "$fastqc_output" -o "$multiqc_output"