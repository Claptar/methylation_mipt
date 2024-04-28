#!/bin/bash

# Load singularity
module load singularity/3.6.1

# Methylatioin calling step
for infile in bismarkOutput/*/*deduplicated.bam
do
    base=$(basename -a -s _trimmed_1_bismark_bt2_pe.deduplicated.bam ${infile})
    echo ${base}
    singularity run singularityImages/bismark.sif bismark_methylation_extractor -p --output_dir bismarkOutput/${base} \
    --no_overlap --gzip --bedGraph --parallel 10 $infile
done