#!/bin/bash
#SBATCH --job-name=your_trait_loci_lister
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=2
#SBATCH --output=example-data_test/logs/p09_locus_breaker_%a.log
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-5 ### depending on the number of you traits

eval "$(conda shell.bash hook)"
conda activate coloc_pipe_env ### include path to env if necessary

### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript --vanilla prj_008_multi_coloc_dev/scripts/locus_breaker_wrap.R \
--path "./prj_008_multi_coloc_dev/example-data" \
--pref "" \
--suf ".1_10546866_11576788.tsv" \
--out "./example-data_test/p09_locus_breaker" \
--sig_pval 5e-8 \
--limit 1e-5 \
--hole 250000 \
--chr "CHROM" \
--pos "GENPOS" \
--pvalue "LOG10P"
