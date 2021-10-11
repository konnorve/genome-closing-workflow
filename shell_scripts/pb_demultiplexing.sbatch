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
#SBATCH -o logs/pb_demul/%j_slurm_output.txt
#SBATCH -e logs/pb_demul/%j_slurm_error.txt

demul_out_path="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch2/analysis/pacbio/cromwell_demul/"
barcodes_fasta_path="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch1/raw/raw_barcodes/Sequel_384_barcodes_v1.fasta"
barcodes_xml_path="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch1/raw/raw_barcodes/Sequel_384_barcodes_v1.xml"
raw_reads_path="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch2/raw/m54030R1_210923_175726.subreadset.xml"

dataset create --force --name my_barcodes --type BarcodeSet ${barcodes_xml_path} ${barcodes_fasta_path}

pbcromwell run pb_demux_subreads -e ${raw_reads_path} -e ${barcodes_xml_path} --output-dir ${demul_out_path} --overwrite
