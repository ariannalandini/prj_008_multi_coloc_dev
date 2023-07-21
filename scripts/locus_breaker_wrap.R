#!/usr/bin/env Rscript
library(data.table)
library(optparse)
library(dplyr)

### Function
source("prj_008_multi_coloc_dev/scripts/loci_identification_funs.R")

### Array job
num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

### Arguments
option_list <- list(
  make_option("--path", type="character", default=NULL,
              help="Path to directory(ies) having all munged GWAS summary statistics", metavar="character"),
  make_option("--pref", type="character", default="", 
              help="Prefix to the trait in the file name (will be removed)", metavar="character"),  
  make_option("--suf", type="character", default=".regenie.gz", 
              help="Suffix to the trait in the file name (will be removed)", metavar="character"),  
  make_option("--out", type="character", default=NULL, 
              help="Path and prefix of output files", metavar="character"),
  make_option("--sig_pval", type="numeric", default=5e-08, 
              help="Significant p-value threshold for top hits", metavar="numeric"),
  make_option("--limit", type="numeric", default=1e-05, 
              help="P-value threshold for loci borders", metavar="numeric"),  
  make_option("--hole", type="integer", default=250000,
              help="Minimum pair-base distance between SNPs in different loci", metavar="integer"),
  make_option("--chr", type="character", default="CHROM",
              help="Name of chromosome column in GWAS summary statistics", metavar="integer"),
  make_option("--pos", type="character", default="GENPOS", 
              help="Name of genetic position column in GWAS summary statistics", metavar="character"),  
  make_option("--pvalue", type="character", default="LOG10P", 
              help="Name of p-value of effect column in GWAS summary statistics", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

########### to delete
#opt$path="/group/pirastu/prj_004_variant2function/gwas_topmed_rap/sum_stats,gwas_hrc_local/ukbb_astle_2016_blood_traits_normalised/results"
#opt$pref="ukbb_topmed_"
#opt$suf="_normalised.tsv.gz"
#opt$out="output_dir1"
###########

## Throw error message - munged GWAS summary statistics file file MUST be provided!
if(is.null(opt$path)){
  print_help(opt_parser)
  stop("Please specify the path of your GWAS summary statistics in --path option", call.=FALSE)
}

### Split list of input paths (divided)
opt$path_list <- unlist(strsplit(opt$path, ","))


### Select one from the list of all available GWAS at the listed path (based on array job) - getting only files ending in ".gz" and ".tsv". List to expand??
gwas <- c(
  list.files(path=opt$path_list, pattern="*.gz$", full.names = T),
  list.files(path=opt$path_list, pattern=".tsv$", full.names = T)
)[[num]]

### Create output folder
if(is.null(opt$out)){ opt$out <- paste0(gsub("^(.*)/.*/.*$", "\\1", gwas), "/loci_definition") } 
system(paste0("mkdir -p ", opt$out))  


### Load-in munged sum stat
sum_stat <- fread(gwas, data.table=F)

### Check if p-value column provided is log10 transformed. If yes, compute original p-value
if (all(sum_stat[[opt$pvalue]] >= 0 & sum_stat[[opt$pvalue]] <= 1)) {
  opt$pvalue <- opt$pvalue
} else {
  sum_stat <- sum_stat %>%
    mutate(P=10^(-sum_stat[[opt$pvalue]]))
  opt$pvalue="P"
}
  
### Loci identification
trait.res <- locus.breaker(sum_stat,
              p.sig=opt$sig_pval,
              p.limit=opt$limit,
              hole.size=opt$hole,
              p.label=opt$pvalue,
              chr.label=opt$chr,
              pos.label=opt$pos)

### Extract trait name. Discard option(s) where "/" is in the name
trait_name <- sapply(opt$path_list, function(x){gsub(paste0(x, "/", opt$pref, "(.*)", opt$suf), "\\1", gwas)}, USE.NAMES=FALSE)
trait_name <- trait_name[which(!grepl("/", trait_name))]

### Add trait name and GWAS sum stat path to the loci table (needed later for coloc). Save
trait.res <- trait.res %>% mutate(trait=trait_name, path=gwas)
cat(paste0("\n", nrow(trait.res), " loci identified for ", trait_name, "\n"))
fwrite(trait.res, paste0(opt$out, "/", trait_name, "_loci.tsv"), sep="\t", quote=F, na=NA)

cat("\n**** DONE!! ****\n")
