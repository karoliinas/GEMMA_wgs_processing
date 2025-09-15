#!/usr/bin/env bash
#
#
# Run somalier on a multisample vcf to infer sample relatedness
# Ped -file format (tab separated) for each sample:
# fam   kid dad mom sex pheno
#
# somalier instructions: https://github.com/brentp/somalier
# somalier sites file downloaded from: https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
# somalier labels downloaded from: https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
# somalier_1kg downloaded from: https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz
#
# Usage: ./somalier.sh multisample.vcf.gz
#
export vcf=$1
export sites="/mnt/gemma/bin/resources/somalier/sites.hg38.vcf.gz"
export labels="/mnt/gemma/bin/resources/somalier/ancestry-labels-1kg.tsv"
export somalier_1kg="/mnt/gemma/bin/resources/somalier/1kg-somalier/*.somalier"
export ref="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export ped="/mnt/data/somalier/gemma.fam"

somalier extract -d extracted/ --sites $sites -f $ref $vcf
somalier relate --ped $ped extracted/*.somalier
