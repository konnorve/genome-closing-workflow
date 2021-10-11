
FLYE_ASSEMBLY_DIR="/nobackup1c/users/kve/2021_PacBioGenomeClosing/batch1/analysis/alternative/flye_assembly"
FLYE_LINK_DIR="/nobackup1/kve/2021_PacBioGenomeClosing/batch1/analysis/alternative/assem_metadata_links"

for barcode in "bc1012" "bc1016" "bc1017" "bc1018" "bc1019" "bc1020" "bc1021" "bc1022"
# barcode="bc${SLURM_ARRAY_TASK_ID}"
do

    echo ${barcode}

    flye_assembly=${FLYE_ASSEMBLY_DIR}/${barcode}/"assembly_info.txt"
    flye_link=${FLYE_LINK_DIR}/assembly_info.${barcode}.tsv

    ln ${flye_assembly} ${flye_link}

done