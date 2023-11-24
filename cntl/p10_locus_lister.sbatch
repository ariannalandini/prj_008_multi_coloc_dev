#!/bin/bash
#SBATCH --job-name=locus_lister
#SBATCH --output=example-data_test/logs/p10_locus_lister.log
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00

eval "$(conda shell.bash hook)"
conda activate coloc_pipe_env ### include path to env if necessary

### The --vanilla option tells Rscript to run without saving or restoring anything in the process. This keeps things nice and clean
Rscript prj_008_multi_coloc_dev/scripts/locus_lister_wrap.R \
--loci_path "./example-data_test/p09_locus_breaker" \
--out "./example-data_test/p10_locus_lister/example-data_all_loci.tsv"