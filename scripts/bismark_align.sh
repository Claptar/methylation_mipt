#!/bin/bash
#PBS -l walltime=10:30:00
#PBS -l nodes=1:ppn=8

### Change directory
#cd /gss/home/a.binkevich.ext/projects/compbio/methylation
SAMPLE=$1

### Make output dir
mkdir -p bismarkOutput/${SAMPLE}

### Load singularity
module load singularity/3.6.1 

### Alignment
singularity run singularityImages/bismark.sif bismark --parallel 4 -p 4  --gzip  references  -1 fastq_trimmed2/${SAMPLE}_trimmed_1.fastq.gz  -2 fastq_trimmed2/${SAMPLE}_trimmed_2.fastq.gz --output_dir bismarkOutput/${SAMPLE} &> bismarkOutput/${SAMPLE}/bismark.log 