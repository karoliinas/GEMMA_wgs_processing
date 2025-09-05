#!/usr/bin/env bash
#
# GEMMA germline variant calling on bwa-mem aligned bam -files using deepvariant-1.6.1
# Change the number of cpus to suit the available resources
#
# Outputs all logs and variant statistics to dir named after bam -file prefix (sampleid)
#
#
#
export cpus=12
export genome_ref="/mnt/gemma/bin/resources/homo_sapiens/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export docker_mnt="/mnt:/mnt"
export regions="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
export run_file="run_DV.cmd"
export wd=$(pwd)
for file in *.bam;
do
    export sample=$(echo $file | awk -F'.' '{print $1}')
    echo "Call variants for $sample"
    export reads=$wd"/"$file
    export vcf=$wd"/"$sample".vcf.gz"
    export gvcf=$wd"/"$sample".g.vcf.gz"
    export intermediate="${wd}/${sample}_intermediate_files"
    export log="${wd}/${sample}_logs"
    echo $reads $vcf $gvcf
    echo $log $intermediate
    echo "# $sample" >> $run_file
    echo "docker run -v $docker_mnt google/deepvariant /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$genome_ref --reads=$reads --sample_name=$sample --output_vcf=$vcf --output_gvcf=$gvcf --intermediate_results_dir=$intermediate --regions '$regions' --num_shards=$cpus --logging_dir=$log" >> $run_file
    echo "" >> $run_file
done


./$run_file
