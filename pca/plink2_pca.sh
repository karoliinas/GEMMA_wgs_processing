#!/usr/bin/env bash
# 
# Convert GEMMA data to plink2 format and run PCA on unrelated samples
#
# ./plink2_pca.sh input_variants.vcf.gz family_pedigree.ped
#
export vcf=$1
export ped=$2
export bfile="${vcf/.vcf.gz/_plink2}'
# Convert VCF to plink format
/mnt/gemma/apps/plink2 \
  --vcf $vcf \
  --make-bed \
  --fam $ped \
  --out $bfile \
  --vcf-half-call m \
  --keep-allele-order \
  --double-id \
  --allow-extra-chr \
  --vcf-idspace-to _ \
  --set-all-var-ids @:# \
  --rm-dup force-first \
  --threads 16 \
  --split-par hg38

export unrelated=$bfile"_unrelated"
/mnt/gemma/apps/plink2 --bfile $bfile --king-cutoff 0.125 --out unrelated


/mnt/gemma/apps/plink2 \
  --bfile $bfile \
  --keep ${bfile}.king.cutoff.in.id \
  --hwe 1e-6 midp \
  --geno 0.1 \
  --make-bed \
  --out unrelated_filt


/mnt/gemma/apps/plink2 \
  --bfile unrelated_filt \
  --indep-pairwise 50 5 0.2 \
  --out ldpruned


#PCA on unrelated samples to generate a population model -> project related samples to it
/mnt/gemma/apps/plink2 \
  --bfile unrelated_filt \
  --keep unrelated.king.cutoff.in.id \
  --extract ldpruned.prune.in \
  --pca allele-wts vcols=chrom,ref,alt \
  --out pca_unrelated

/mnt/gemma/apps/plink2 \
  --bfile $bfile \
  --extract ldpruned.prune.in \
  --score unrelated.eigenvec.var 2 4 header-read no-mean-imputation variance-standardize  \
  --score-col-nums 5-14 \
  --out pca_all_projected
