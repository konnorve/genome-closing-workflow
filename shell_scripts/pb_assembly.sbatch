#!/usr/bin/env bash
# 
#SBATCH --job-name=pb_assem
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 10
#
#SBATCH --comment="pacbio assembly"
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
#SBATCH -o logs/pb_assem/%j_%a_slurm_output.txt
#SBATCH -e logs/pb_assem/%j_%a_slurm_error.txt

demul_dir="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch2/analysis/pacbio/cromwell_demul/outputs"
assem_dir="/nobackup1/kve/2021_HIFI_Genome_Closing_Project/batch2/analysis/pacbio/pb_assembly"

# batch 1: "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
# batch 2: "bc1015" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
barcode="bc${SLURM_ARRAY_TASK_ID}"

echo ${barcode}

subread_set=${demul_dir}/"demultiplex.${barcode}--${barcode}.subreadset.xml"
out_dir=${assem_dir}/${barcode}

mkdir -p ${out_dir}

echo ${subread_set}
echo ${out_dir}

pbcromwell run pb_assembly_microbial -e ${subread_set} --task-option microasm_genome_size=2000000 --task-option microasm_coverage=30 --output-dir ${out_dir} --overwrite