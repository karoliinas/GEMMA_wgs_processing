For fastq in *1.fq.gz;do
  R1_IN=$fastq
  R2_IN="${fastq/1.fq.gz/2.fq.gz}
  R1_OUT="${R1_IN/1.fq.gz/1.trimmed.fq.gz}"
  R2_OUT="${R2_IN/2.fq.gz/2.trimmed.fq.gz}"

  # Run fastp (trim adapters, poly g and overrepresented sequences)
  fastp \
    -i $R1_IN -I $R2_IN \
    -o $R1_OUT -O $R2_OUT \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --overrepresentation_analysis
