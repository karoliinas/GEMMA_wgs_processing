#!/usr/bin/env bash
#
#
# Run bcftools stats on a set of variants for a multisample vcf
# after joint genotyping with GLnexus
#
#
# usage: ./stats.sh joint_vcf.vcf.gz
#
export genomeref="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna gemma239_norm.vcf.gz"
export vcf=$1
export stats="${vcf/.vcf.gz/.stats}"

bcftools stats -F $genomeref $vcf > $stats
bcftools plot-vcfstats -s $stats

