#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 48:00:00
#SBATCH -J I2
#SBATCH -o OUT/I2.out
#SBATCH -e ERR/I2.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
spaceranger count --id=I2 \
    --transcriptome=ssc97Spac \
    --fastqs=FastQ \
    --sample=I2-STB \
    --image=JsonTIF/IPP_STB.tif \
    --slide=V10M24-003 \
    --area=B1 \
    --loupe-alignment=JsonTIF/IPP_STB.json \
    --localcores=38

