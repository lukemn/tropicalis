#!/bin/bash

module purge
source ~/py3.6.3/bin/activate

module load bcftools/intel/1.9
module load samtools/intel/1.9
module load bedtools/intel/2.27.1
#module load bwa/intel/0.7.17
module load r/intel/3.4.2
module load htslib/intel/1.4.1
module load gatk/4.1.3.0
module load picard/2.17.11

export PATH=$PATH:~/bin/bwa-mem2-2.0pre2_x64-linux/

alias sk="snakemake --cores $SLURM_CPUS_PER_TASK"


