#!/bin/bash
#SBATCH -p tsl-long                             # queue
#SBATCH --wckey=talbot_core                    # project_code
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 26                                   # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH -o slurm.%j.out                         # slurm job output
#SBATCH -e slurm.%j.err                         # slurm job error
#SBATCH --mail-type=begin,end,fail              # notifications for job start, end & fail
#SBATCH --mail-user=neha.sahu@tsl.ac.uk         # send-to address

source trimmomatic-0.36

# Input directory containing FASTQ files
input_dir1="/tsl/data/reads/ntalbot/timecourse_rnaseq_bip1_mst12/"
#input_dir2="/tsl/data/reads/ntalbot/guy11_rnaseq_appressorium"

# Output directory for trimmed files
output_dir="trimmed"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Run script over each sample name (assuming the sample names are unique up to the last "_")
find "$input_dir1" -name "*_1.fq.gz" | while read -r file; do
    # Extract the sample name without the suffix "_1"
    sample_name="${file%_1.fq.gz}"

    # Input file paths for both pairs
    input_file_1="${sample_name}_1.fq.gz"
    input_file_2="${sample_name}_2.fq.gz"

    # Output file paths
    output_file_paired_1="$output_dir/$(basename "${sample_name}")_paired_1.fastq.gz"
    output_file_unpaired_1="$output_dir/$(basename "${sample_name}")_unpaired_1.fastq.gz"
    output_file_paired_2="$output_dir/$(basename "${sample_name}")_paired_2.fastq.gz"
    output_file_unpaired_2="$output_dir/$(basename "${sample_name}")_unpaired_2.fastq.gz"
    log_file="$output_dir/$(basename "${sample_name}")_trimmomatic.log"

    # Run Trimmomatic for paired-end reads
    trimmomatic PE -phred33 "$input_file_1" "$input_file_2" \
        "$output_file_paired_1" "$output_file_unpaired_1" \
        "$output_file_paired_2" "$output_file_unpaired_2" \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:25 \
        ILLUMINACLIP:/tsl/software/testing/trimmomatic/0.36/x86_64/share/trimmomatic/adapters/ilmn_adapters.fa:2:40:15 \
        > "$log_file" 2>&1  # output log file
done