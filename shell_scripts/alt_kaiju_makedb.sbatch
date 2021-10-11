#!/usr/bin/env bash
# 
#SBATCH --job-name=kaiju_makedb
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
#SBATCH -o /home/kve/scripts/pb_pro_assembly/genome-closing-workflow/shell_scripts/logs/kaiju_makedb/%j_%a_slurm_output.txt
#SBATCH -e /home/kve/scripts/pb_pro_assembly/genome-closing-workflow/shell_scripts/logs/kaiju_makedb/%j_%a_slurm_error.txt

KAIJU_DB_DIR="/pool001/kve/kaiju_db"

cd ${KAIJU_DB_DIR}

kaiju-makedb -t 10 -s refseq
