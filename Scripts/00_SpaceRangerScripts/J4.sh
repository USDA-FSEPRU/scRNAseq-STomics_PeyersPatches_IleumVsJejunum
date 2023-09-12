#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J J4
#SBATCH -o OUT/J4.out
#SBATCH -e ERR/J4.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=J4 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=J4-STD \
    --image=JsonTIF/JPP_STD.tif \
    --slide=V10M24-004 \
    --area=D1 \
    --loupe-alignment=JsonTIF/JPP_STD.json \
    --localcores=38

