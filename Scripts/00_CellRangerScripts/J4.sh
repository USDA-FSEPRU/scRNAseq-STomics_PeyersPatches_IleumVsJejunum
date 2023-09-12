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
cellranger count --id=J4 \
                   --transcriptome=/project/fsepru/PIGscRNAseq1/ssc97Cell \
                   --fastqs=FastQ \
                   --sample=J-4-1,J-4-2,J-4-3,J-4-4 \
                   --localcores=38
