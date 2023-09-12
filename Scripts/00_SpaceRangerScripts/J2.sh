#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J J2
#SBATCH -o OUT/J2.out
#SBATCH -e ERR/J2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=J2 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=J2-STC \
    --image=JsonTIF/JPP_STC.tif \
    --slide=V10M24-004 \
    --area=C1 \
    --loupe-alignment=JsonTIF/JPP_STC.json \
    --localcores=38

