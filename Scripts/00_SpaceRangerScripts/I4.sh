#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J I4
#SBATCH -o OUT/I4.out
#SBATCH -e ERR/I4.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=I4 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=I4-STD \
    --image=JsonTIF/IPP_STD.tif \
    --slide=V10M24-003 \
    --area=D1 \
    --loupe-alignment=JsonTIF/IPP_STD.json \
    --localcores=38

