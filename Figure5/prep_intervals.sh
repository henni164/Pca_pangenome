#!/bin/bash -l
#SBATCH --job-name=bigplot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=50GB
#SBATCH --mail-user=eva.henningsen@csiro.au
#SBATCH --mail-type=FAIL
#SBATCH --account=OD-221325

#to submit: sbatch -a 0-31 --export INPUTLIST="/scratch3/hen294/cactus_variants_pangenome_plot/haps.txt" prep_intervals.sh

module load R

if [ -s "$INPUTLIST" ]
then
    SAMPLES=( $(cut -f 1 $INPUTLIST) );

    if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        i=$SLURM_ARRAY_TASK_ID

	Rscript prep_intervals.R -i ${SAMPLES[i]}

    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
    fi
else
    echo "Error: Missing inputs file list as --export env INPUTLIST or file empty"
fi
