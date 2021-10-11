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
#SBATCH --array=1012,1016,1017,1018,1019,1020,1021,1022
#
#SBATCH -o logs/alt_assem_canu/%j_%a_slurm_output.txt
#SBATCH -e logs/alt_assem_canu/%j_%a_slurm_error.txt

lima_dir="/nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/alternative/lima_demultiplexed"
canu_dir="/nobackup1c/users/kve/2021_PacBioGenomeClosing/batch1/analysis/alternative/canu_assembly"

# "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
barcode="bc${SLURM_ARRAY_TASK_ID}"

echo ${barcode}

fastq_path=${lima_dir}/"lima_demultiplexed.${barcode}.fastq"
canu_bc_dir=${canu_dir}/${barcode}

mkdir ${canu_bc_dir}

echo ${fastq_path}
echo ${canu_bc_dir}

canu -p ${barcode} -d ${canu_bc_dir} \
    genomeSize=2.5m \
    -pacbio ${fastq_path} 

