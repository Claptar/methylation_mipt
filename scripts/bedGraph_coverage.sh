#!/bin/bash

for file in bismarkOutput/*/*.deduplicated.bam
do
    base=$(basename -a -s _trimmed_1_bismark_bt2_pe.deduplicated.bam ${file})
    singularity run singularityImages/samtools.sif samtools sort -@16 $file -o bismarkOutput/${base}/${base}.deduplicated.sorted.bam
    singularity run singularityImages/bedtools.sif bedtools genomecov \
    -ibam bismarkOutput/${base}/${base}.deduplicated.sorted.bam -bg > bismarkOutput/${base}/${base}.coverage.bedGraph
done