#!/bin/bash

# Load singularity
module load singularity/3.6.1

# Dedupliaction step
for infile in bismarkOutput/*/*_trimmed_1_bismark_bt2_pe.bam
do
    base=$(basename -a -s _trimmed_1_bismark_bt2_pe.bam ${infile})
    echo ${base}
    singularity run singularityImages/bismark.sif deduplicate_bismark -p --output_dir bismarkOutput/${base} $infile
done