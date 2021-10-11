#!/usr/bin/env bash
# 
#SBATCH --job-name=alt_kaiju
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
#SBATCH --comment="Kaiju on Flye Assembly"
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
#SBATCH -o logs/alt_kaiju_classification/%j_%a_slurm_output.txt
#SBATCH -e logs/alt_kaiju_classification/%j_%a_slurm_error.txt

ALT_ASSEMBLY_DIR="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/alt_assembly_links"
KAIJU_CLASSIFICATION_DIR="/nobackup1/kve/2021_PacBioGenomeClosing/batch2/analysis/alternative/kaiju_classification"
KAIJU_DB_DIR="/pool001/kve/kaiju_db"

# batch 1: "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
# batch 2: "bc1015" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
barcode="bc${SLURM_ARRAY_TASK_ID}"

echo ${barcode}

assembly_path=${ALT_ASSEMBLY_DIR}/"final_assembly.${barcode}.fasta"
kaiju_bc_result_txt=${KAIJU_CLASSIFICATION_DIR}/"kaiju_classification_output.${barcode}.txt"
kaiju_bc_result_table=${KAIJU_CLASSIFICATION_DIR}/"kaiju_classification_output.${barcode}.tsv"
kaiju_bc_result_names=${KAIJU_CLASSIFICATION_DIR}/"kaiju_classification_output_named.${barcode}.tsv"
kaiju_index_file=${KAIJU_DB_DIR}/"refseq/kaiju_db_refseq.fmi"
nodes_dmp=${KAIJU_DB_DIR}/"nodes.dmp"
names_dmp=${KAIJU_DB_DIR}/"names.dmp"

echo ${assembly_path}
echo ${kaiju_bc_dir}

# kaiju -v -t ${nodes_dmp} -f ${kaiju_index_file} -i ${assembly_path} -o ${kaiju_bc_result_txt} -z ${SLURM_CPUS_ON_NODE}

# kaiju2table -t ${nodes_dmp} -n ${names_dmp} -r genus -o ${kaiju_bc_result_table} ${kaiju_bc_result_txt}

kaiju-addTaxonNames -t ${nodes_dmp} -n ${names_dmp} -i ${kaiju_bc_result_txt} -o ${kaiju_bc_result_names}