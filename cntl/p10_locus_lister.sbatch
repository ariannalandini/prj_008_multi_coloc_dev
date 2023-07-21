#!/bin/bash
#SBATCH --job-name=locus_lister
#SBATCH --output=logs/p10_locus_lister.log
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00


# Loading modules
module load R/4.1.0

### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript prj_008_multi_coloc_dev/scripts/locus_lister_wrap.R \
--loci_path "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/mh_and_loci" \
--snp "ID" \
--chr "CHROM" \
--pos "GENPOS" \
--a1 "ALLELE1" \
--a0 "ALLELE0" \
--effect "BETA" \
--se "SE" \
--pval "P" \
--out "gwas_topmed_rap/mh_and_loci/ukbb_topmed_all_loci.tsv"