#!/usr/bin/env bash
# 
#SBATCH --job-name=assem_smk
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
#SBATCH --comment="Lima demultiplexing"
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
#SBATCH -o logs/alt_lima/%j_slurm_output.txt
#SBATCH -e logs/alt_lima/%j_slurm_error.txt

lima_out_path="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/lima_demultiplexed/lima_demultiplexed.bam"
barcodes_path="/nobackup1/kve/2021_PacBioGenomeClosing/batch1/raw/raw_barcodes/Sequel_384_barcodes_v1.fasta"
raw_reads_path="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/raw/m54030R1_210923_175726.subreadset.xml"


lima --same --peek-guess --split-named --num-threads 10 ${raw_reads_path} ${barcodes_path} ${lima_out_path}