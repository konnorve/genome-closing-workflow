#!/usr/bin/env bash
# 
#SBATCH --job-name=pb_gorg
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
#SBATCH --comment="PacBio GORG Run"
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
# SBATCH --array=1012,1016,1017,1018,1019,1020,1021,1022
#
#SBATCH -o logs/pb_gorg/%j_%a_slurm_output.txt
#SBATCH -e logs/pb_gorg/%j_%a_slurm_error.txt

PB_ASSEMBLY_DIR="/nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/pacbio/assembly"
GORG_DIR="/nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/pacbio/pb_GORG"

for barcode in "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
# barcode="bc${SLURM_ARRAY_TASK_ID}"
do
    echo ${barcode}

    assembly_path=${PB_ASSEMBLY_DIR}/${barcode}/"outputs/assembly.rotated.polished.renamed.fsa"
    gorg_bc_dir=${GORG_DIR}/${barcode}

    mkdir ${gorg_bc_dir}

    echo ${assembly_path}
    echo ${gorg_bc_dir}

    nextflow run BigelowLab/gorg-classifier -profile conda --seqs ${assembly_path} --outdir ${gorg_bc_dir} --cpus ${SLURM_CPUS_ON_NODE}
done