#!/bin/bash
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=1

### Change directory
#cd /gss/home/a.binkevich.ext/projects/compbio/methylation
SAMPLE=$1

### Load singularity
module load singularity/3.6.1 

### Trimming
singularity run singularityImages/fastp.sif fastp -i fastq_files/${SAMPLE}_1.fastq.gz -I fastq_files/${SAMPLE}_2.fastq.gz \
-o fastq_trimmed2/${SAMPLE}_trimmed_1.fastq.gz -O fastq_trimmed2/${SAMPLE}_trimmed_2.fastq.gz \
-h qc/${SAMPLE}_fastp.html -j qc/${SAMPLE}_fastp.json --average_qual 20 --length_required 20 -w 16