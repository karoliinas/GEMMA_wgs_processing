#!/bin/bash
# ============================================================
# Filter VEP-annotated multisample VCF
#   - Keep protein-coding + promoter variants
#   - Restrict to SFARI genes
#   - Apply AF filter (rare < 0.05)
#   - Genotype quality (GQ > 20)
#   - Call rate filter (≥ 80%)
#   - Keep only polymorphic sites
# ============================================================

VCF=$1         # VEP-annotated multisample VCF
GENE_LIST="sfari25.txt"    # One gene symbol per line
AF_FIELD="gnomADg_AF"           # Change to the AF field in your VCF (e.g. AF, gnomAD_AF, MAX_AF)
GQ_MIN=20                      # Minimum genotype quality
CALLRATE=0.8                   # Minimum call rate (proportion of samples with non-missing calls)
OUT="${VCF/vcf.gz/filt.vcf.gz}"     # Output file
CONSEQ="missense_variant|stop_gained|stop_lost|synonymous_variant|frameshift_variant|inframe_insertion|inframe_deletion|promoter_variant"


bcftools +split-vep $VCF --columns SYMBOL,Consequence,$AF_FIELD \
  | bcftools filter -i "FMT/GQ>$GQ_MIN" -S . \                         # set low-GQ genotypes missing
  | bcftools view -c 1 \                                               # remove invariant sites
  | bcftools view -i "F_MISSING<$(echo "1-$CALLRATE" | bc -l) && INFO/$AF_FIELD<0.05 && Consequence ~ \"$CONSEQ\" && SYMBOL=@$GENE_LIST" \
  -Oz -o $OUT

# Index the final VCF
bcftools index -t $OUT

echo "✅ Filtering complete. Output written to: $OUT"

