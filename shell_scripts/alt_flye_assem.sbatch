#!/usr/bin/env bash
# 
#SBATCH --job-name=assem_alt
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
#SBATCH --comment="fly assembly"
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
#SBATCH --array=1015,1016,1017,1018,1019,1020,1021,1022
#
#SBATCH -o logs/alt_assem_flye/%j_%a_slurm_output.txt
#SBATCH -e logs/alt_assem_flye/%j_%a_slurm_error.txt

lima_dir="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/lima_demultiplexed"
fly_dir="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/flye_assembly"

# batch 1: "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
# batch 2: "bc1015" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
barcode="bc${SLURM_ARRAY_TASK_ID}"

echo ${barcode}

fastq_path=${lima_dir}/"lima_demultiplexed.${barcode}.fastq"
fly_bc_dir=${fly_dir}/${barcode}

mkdir -p ${fly_bc_dir}

echo ${fastq_path}
echo ${fly_bc_dir}

flye --pacbio-raw ${fastq_path} --out-dir ${fly_bc_dir}