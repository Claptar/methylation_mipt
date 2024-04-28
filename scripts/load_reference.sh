#!/bin/bash

# скачаем сборку 17-ой хромосомы
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.17.fa.gz
gzip -d Mus_musculus.GRCm39.dna.chromosome.17.fa.gz
mv Mus_musculus.GRCm39.dna.chromosome.17.fa references/chr17.fa

# скачаем геномную аннотацию с ENSEMBL
wget https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz
gzip -d Mus_musculus.GRCm39.111.gtf.gz
mv Mus_musculus.GRCm39.111.gtf references/genome.gtf