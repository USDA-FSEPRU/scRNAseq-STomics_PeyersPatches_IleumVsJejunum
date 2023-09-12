#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J J3
#SBATCH -o OUT/J3.out
#SBATCH -e ERR/J3.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=J3 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=J3-STB \
    --image=JsonTIF/JPP_STB.tif \
    --slide=V10M24-004 \
    --area=B1 \
    --loupe-alignment=JsonTIF/JPP_STB.json \
    --localcores=38

