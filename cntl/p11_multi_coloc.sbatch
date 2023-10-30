#!/bin/bash
#SBATCH --job-name=p11_locus_%a.log
#SBATCH --output=coloc/multi_coloc/logs/p11_multi_coloc_locus_%a.log ### fix this - output folder should be specified in the Rscript?
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=20:00:00
#SBATCH --array=1-1328%100

# Loading modules
module load R/4.1.0
module load bcftools/1.14-python-3.9.10 ### required for tabix


### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript --vanilla prj_008_multi_coloc_dev/scripts/multi_coloc_wrap.R \
--input "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/mh_and_loci/ukbb_topmed_all_loci.tsv" \
--output "coloc/multi_coloc" \
--chr "CHROM" \
--pos "GENPOS" \
--rsid "ID" \
--a1 "ALLELE1" \
--a0 "ALLELE0" \
--freq "A1FREQ" \
--n "N" \
--effect "BETA" \
--se "SE" \
--pvalue "P" \
--grch 38 \
--maf 0.0001 \
--bfile path/to/my/custom_ld/name_of_plink_file \
--save_inter_files "TRUE"
