#!/bin/bash
#PBS -N quality_control
#PBS -l walltime=10:30:00
#PBS -l nodes=1:ppn=8

### Change directory
cd /gss/home/a.binkevich.ext/projects/compbio/methylation

### Load singularity
module load singularity/3.6.1 

### Run QC
singularity run singularityImages/fastqc.sif fastqc -o qc -t 8 fastq_files/*.fastq.gz
singularity run singularityImages/multiqc.sif multiqc qc --outdir qc --filename multiqc.html --title "Read quality" --force