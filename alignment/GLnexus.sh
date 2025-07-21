#!/bin/bash
#
#
#
# 
export input_dir="/mnt/data/variants"
export bed="/data/references/GCA_000001405.15_GRCh38_no_alt_analysis_set_mainChr.bed" #regions to process
export output1="/data/variants/GLnexus/GEMMA.bcf"
export output2="/data/variants/GLnexus/GEMMA.vcf.gz"
export output3="/data/variants/GLnexus/GEMMA_norm.vcf.gz"
export log="/data/variants/GLnexus/GEMMA.log"

docker run --security-opt seccomp=unconfined -v /mnt:/mnt ghcr.io/dnanexus-rnd/glnexus:v1.4.1 glnexus_cli --config DeepVariant \
  --dir /mnt/data/variants/GLnexus/db/GLnexus.DB --bed $bed --trim-uncalled-alleles $input_dir/*.g.vcf.gz 1> $output 2> $log
bcftools view $output -Oz -o $output2
tabix $output2
bcftools norm -m - $output2 -Oz -o $output3
tabix $output3
