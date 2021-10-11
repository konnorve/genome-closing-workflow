#!/usr/bin/env bash
# 
#SBATCH --job-name=bam2fq
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 12:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 10
#
#SBATCH --comment="bam2fastq conversion"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
# Request nodes:
#SBATCH --nodes=1
#
#SBATCH --mem 250000
#
#SBATCH -o logs/alt_bam2fastq/%j_slurm_output.txt
#SBATCH -e logs/alt_bam2fastq/%j_slurm_error.txt

for barcode in "bc1015" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
do
echo ${barcode}
lima_dir="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/lima_demultiplexed"
bam_path=${lima_dir}/"lima_demultiplexed.${barcode}--${barcode}.bam"
fastq_path=${lima_dir}/"lima_demultiplexed.${barcode}.fastq"

echo ${bam_path}
echo ${fastq_path}

samtools fastq ${bam_path} > ${fastq_path}
done
