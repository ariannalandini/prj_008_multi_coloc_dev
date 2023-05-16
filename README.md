# prj_008_multi_coloc_dev

 - Goal of the pipeline
 - Link to paper(s) of interest


## Description
Performs:

1) Identification of trait-specific loci 
2) Identification of pan loci, grouping multiple overlapping loci across traits 
3) Colocalisation with decomposition of association signal 
    3.1) Identification, for each trait at each locus, of indipendent association signals 
    3.2) Leave-one-out conditioning of indipendent SNPs to clear the association signals 
    3.3) Coloc between all possible pairs of indipendent SNP-trait 


## Installation
- Should we create a specific environment with all R packages used/COJO ?


## Usage

1) Create trait-specific loci table
To identify the associated loci specific to each trait in a list, you should run the cntl/p09_locus_breaker.sbatch script, providing:\
    --path: path to directory where GWAS summary statistics for all your traits of interest are located.\
Note this is the only essential option, if not provided the script will throw an error.\
    --pref: prefix to the trait name of GWAS summary statistics\
    --suf: suffix to the trait name of GWAS summary statistics\
Prefix and suffix are necessary so that only the trait name (e.g. "height", "bmi", etc.) is extracted.\
    --out: path to directory where output tables will be stored\
While all the listed options are quite relevant, they are not strictly essential: for example, if your GWAS summary statistics have been created using the [Regenie pipeline](https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie) and you are happy with having a "loci_identification" folder created in the parent (??) of your input folder, then you are fine with defaul arguments.\
    --chr: Name of chromosome column in GWAS summary statistics (default: CHROM)\
    --pos: Name of genetic position column in GWAS summary statistics (default: GENPOS)\
    --pvalue: Name of p-value of effect column in GWAS summary statistics (default: P)\
    --sig_pval: Significant p-value threshold for top hits (default:5e-08)\
    --limit: P-value threshold for loci borders (default: 1e-05)\
    --hole: Minimum pair-base distance between SNPs in different loci (default: 250000)\


2) Create pan-loci table

3) Run coloc


## To do
- Patch coloc.abf() function to account for sample overlap
- Slim down output files?
- Include possibility to perform coloc on a specified custom region
- Include TileDB format as input


## Authors and acknowledgment
Nicola Pirastu (nicola.pirastu@fht.org)\
Arianna Landini (arianna.landini@external.fht.org)\
Sodbo Sharapov (sodbo.sharapov@fht.org)\
