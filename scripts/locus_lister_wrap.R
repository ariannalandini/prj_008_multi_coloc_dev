#### Pietro's script to identify overlapping loci across multiple traits and collapse them in a mega locus

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
source("prj_008_multi_coloc_dev/scripts/loci_identification_funs.R")

### Set arguments
option_list <- list(
  make_option("--loci_path", type="character", default=NULL, 
              help="Path to directory having all precomputed loci", metavar="character"),
  make_option("--out", type="character", default="./panlocus_table.tsv", 
              help="Path and name of output directory and file", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### Run locus_lister!
locus.lister(my_path_loci=opt$loci_path, out=opt$out)
cat("\n**** DONE!! ****\n")