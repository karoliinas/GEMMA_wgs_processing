#!/bin/bash
#
# Annotate variants with vep ->
#     * all available annotations
#     * pick the most common
#     * VEP cache files downloaded from https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
#     * CADD scores downloaded from https://cadd.gs.washington.edu/download

export vcf=$1
export output=${vcf/.vcf.gz/_VEP.vcf}"
export reference="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export cache_dir="/mnt/gemma/bin/resources/vep/"
export cadd_snv="/mnt/gemma/bin/resources/vep/whole_genome_SNVs.tsv.gz"
export cadd_indel="/mnt/gemma/bin/resources/vep/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz"


docker run --security-opt seccomp=unconfined -v /mnt:/mnt ensemblorg/ensembl-vep vep \
  -i $vcf -o $output \
  --fasta $reference --cache --offline --format vcf \
  --vcf --everything --fork 4 --pick --dir_cache $cache_dir \
  --plugin CADD,snv=$cadd_snv,indels=$cadd_indel
