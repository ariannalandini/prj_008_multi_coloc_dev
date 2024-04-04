#!/bin/bash
#SBATCH --job-name=p11_locus_%a.log
#SBATCH --output=example-data_test/logs/p11_multi_coloc_locus_%a.log ### fix this - output folder should be specified in the Rscript?
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=30:00:00
#SBATCH --array=1
#1-1328%100

eval "$(conda shell.bash hook)"
conda activate coloc_pipe_env ### include path to env if necessary

### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript --vanilla prj_008_multi_coloc_dev/scripts/multi_coloc_wrap.R \
--input "./example-data_test/p10_locus_lister/example-data_all_loci.tsv" \
--output "./example-data_test/p11_multi_coloc_locus" \
--chr "CHROM" \
--pos "GENPOS" \
--rsid "ID" \
--a1 "ALLELE1" \
--a0 "ALLELE0" \
--freq "A1FREQ" \
--n "N" \
--effect "BETA" \
--se "SE" \
--pvalue "LOG10P" \
--maf 0.0001 \
--pph4 0.9 \
--save_inter_files "TRUE" \
--tabix_bin "/ssu/gassu/software/htslib-tools/1.14/tabix" \
--plink2_bin "/ssu/gassu/software/plink/2.00_20211217/plink2" \
--gcta_bin "/ssu/gassu/software/GCTA/1.94.0beta/gcta64"
