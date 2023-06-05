#### Pietro's script to identify overlapping loci across multiple traits and collapse them in a mega locus

### Must provide: 
# 1) my_path_gwas - path to GWAS
# 2) traits_pref AND traits_sub - prefix and suffix of GWAS files

## Other arguments should have default options - to revisit once established the tileDB format


### TO DO:
# Add if/else statement to check these conditions are met


library(optparse)

### Function
source("/group/pirastu/prj_004_variant2function/scripts/loci_identification_funs.R")

### Set arguments

option_list <- list(
  make_option("--gwas_path", type="character", default=NULL, 
              help="Path to directory having all GWAS summary statistics", metavar="character"),
  make_option("--pref", type="character", default="", 
              help="Instead of providing the list of traits names, you can provide the file name prefix to the trait (to be removed)", metavar="character"),  
  make_option("--suf", type="character", default=".regenie.gz", 
              help="Instead of providing the list of traits names, you can provide the file name suffix to the trait (to be removed)", metavar="character"),
  make_option("--sig_pval", type="numeric", default=5e-08, 
              help="Significant p-value threshold for top hits", metavar="numeric"),
  make_option("--limit", type="numeric", default=1e-05, 
              help="P-value threshold for loci borders", metavar="numeric"),  
  make_option("--hole", type="integer", default=250000,
              help="Minum pair-base distance between SNPs in different loci", metavar="integer"),
  make_option("--loci_path", type="character", default=NULL, 
              help="Path to directory having all precomputed loci", metavar="character"),
  make_option("--snp", type="character", default="ID", 
              help="Name of rsid column", metavar="character"),  
  make_option("--chr", type="character", default="CHROM",
              help="Name of chromosome column", metavar="integer"),
  make_option("--pos", type="character", default="GENPOS", 
              help="Name of genetic position column", metavar="character"),  
  make_option("--a1", type="character", default="ALLELE1", 
              help="Name of effect allele column", metavar="character"),  
  make_option("--a0", type="character", default="ALLELE0", 
              help="Name of NON effect allele column", metavar="character"),  
  make_option("--effect", type="character", default="BETA", 
              help="Name of effect size column", metavar="character"),  
  make_option("--se", type="character", default="SE", 
              help="Name of standard error of effect column", metavar="character"), 
  make_option("--pval", type="character", default="P", 
              help="Name of p-value of effect column", metavar="character"),
  make_option("--out", type="character", default="./panlocus_table.tsv", 
              help="Path and name of output directory and file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


### Run locus_lister!

locus.lister(
  my_path_gwas=opt$path,
  traits_pref=opt$pref,
  traits_suf=opt$suf,
  p.sig=opt$sig_pval,
  p.limit=opt$limit,
  hole.size=opt$hole,
  my_path_loci=opt$loci_path,
  SNP=opt$snp,
  CHR=opt$chr,
  BP=opt$pos,
  A1=opt$a1,
  A2=opt$a0,
  BETA=opt$effect,
  SE=opt$se,
  P=opt$pval,
  out=opt$out)
cat("\n**** DONE!! ****\n")

