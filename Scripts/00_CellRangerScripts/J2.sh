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
cellranger count --id=J2 \
                   --transcriptome=/project/fsepru/PIGscRNAseq1/ssc97Cell \
                   --fastqs=FastQ \
                   --sample=J-2-1,J-2-2,J-2-3,J-2-4 \
                   --localcores=38
