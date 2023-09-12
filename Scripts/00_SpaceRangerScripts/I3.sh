#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J I3
#SBATCH -o OUT/I3.out
#SBATCH -e ERR/I3.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=I3 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=I3-STC \
    --image=JsonTIF/IPP_STC.tif \
    --slide=V10M24-003 \
    --area=C1 \
    --loupe-alignment=JsonTIF/IPP_STC.json \
    --localcores=38

