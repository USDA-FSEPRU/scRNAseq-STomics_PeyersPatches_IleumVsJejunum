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
cellranger count --id=I3 \
                   --transcriptome=/project/fsepru/PIGscRNAseq1/ssc97Cell \
                   --fastqs=FastQ \
                   --sample=I-3-1,I-3-2,I-3-3,I-3-4 \
                   --localcores=38
