#!/bin/bash
#SBATCH --job-name=blood_loci_lister
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/p09_locus_breaker_%a.log
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --array=5-29

# Loading modules
module load R/4.1.0

### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript --vanilla /group/pirastu/prj_004_variant2function/scripts/locus_breaker_wrap.R \
--path "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/sum_stats" \
--pref "ukbb_topmed_" \
--suf "_normalised.tsv.gz" \
--out "/group/pirastu/prj_004_variant2function/gwas_topmed_rap"\

