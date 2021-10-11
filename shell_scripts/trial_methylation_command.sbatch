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
#SBATCH --partition newnodes
#
# Number of tasks/cores for job
#SBATCH -n 10
#
#SBATCH --comment="Genome Closing with multiple assemblers."
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
#SBATCH -o logs/sbatch/%j_slurm_output.txt
#SBATCH -e logs/sbatch/%j_slurm_error.txt

pbcromwell run pb_basemods -e /nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/pacbio/pb_demultiplexed/outputs/demultiplex.bc1020--bc1020.subreadset.xml -e /nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/pacbio/assembly/bc1020/outputs/reference.dataset.xml --output-dir /nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/pacbio/methylation/bc1020 --overwrite
