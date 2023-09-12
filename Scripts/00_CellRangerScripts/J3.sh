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
cellranger count --id=J3 \
                   --transcriptome=/project/fsepru/PIGscRNAseq1/ssc97Cell \
                   --fastqs=FastQ \
                   --sample=J-3-1,J-3-2,J-3-3,J-3-4 \
                   --localcores=38
