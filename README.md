# prj_008_multi_coloc_dev

 - Goal of the pipeline
 - Link to paper(s) of interest


## Description
Performs:

1) Identification of trait-specific loci
2) Identification of pan loci, grouping multiple overlapping loci across traits
3) Colocalisation with decomposition of association signal\
    3.1) Identification, for each trait at each locus, of indipendent association signals\
    3.2) Leave-one-out conditioning of indipendent SNPs to clear the association signals\
    3.3) Coloc between all possible pairs of indipendent SNP-trait


## Installation
- Should we create a specific environment with all R packages used/COJO ?


## Quick start
- Give minimum example (example data needed?)

Specify input by providing argument to the `--path` option in the `p09_locus_breaker.sbatch`
```
Rscript --vanilla /group/pirastu/prj_004_variant2function/scripts/locus_breaker_wrap.R \
--path "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/sum_stats"
```
Then run specifying the number of traits
```
sbatch cntl/p09_locus_breaker.sbatch --array=1-29
```


```
sbatch cntl/p10_locus_lister.sbatch
```
```
Rscript /group/pirastu/prj_004_variant2function/scripts/locus_lister_wrap.R \
--loci_path "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/mh_and_loci"
```

```
sbatch cntl/p11_multi_coloc.sbatch
```
```
Rscript --vanilla /group/pirastu/prj_004_variant2function/scripts/multi_coloc_wrap.R \
--input "/group/pirastu/prj_004_variant2function/gwas_topmed_rap/mh_and_loci/ukbb_topmed_all_loci.tsv"
```


## Usage

1) Create trait-specific loci table

To identify the associated loci specific to each trait in a list, run the **`cntl/p09_locus_breaker.sbatch`** script, providing:\
    `--path`: path to directory or multiple (comma separated) directories where GWAS summary statistics for all your traits of interest are located.\
Note this is the only essential option, if not provided the script will throw an error and stop.\
    `--pref`: prefix to the trait name of provided GWAS summary statistics\
    `--suf`: suffix to the trait name of provided GWAS summary statistics\
Prefix and suffix are necessary to extract only the trait name (e.g. "height", "bmi", etc.) from your GWAS summary statistics file names. Of course, this assumes that all files have been named following the same logic.\
    `--out`: path to the (single) directory where output loci tables will be stored. \
    `--chr`: Name of chromosome column in GWAS summary statistics (default: CHROM)\
    `--pos`: Name of genetic position column in GWAS summary statistics (default: GENPOS)\
While all the listed options are quite relevant, they are not strictly essential: for example, if your GWAS summary statistics have been created using the [GAU Regenie pipeline](https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie) and you are happy with having a `"loci_identification"` directory created in the parent directory of your input folder, then you are fine with defaul arguments.\
    `--pvalue`: Name of p-value of effect column in GWAS summary statistics (default: P)\
NB: Regenie (and thus the [GAU Regenie pipeline](https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie)) will output only the LOG10 p-value
    `--sig_pval`: Significant p-value threshold for top hits (default:5e-08)\
    `--limit`: P-value threshold for loci borders (default: 1e-05)\
    `--hole`: Minimum pair-base distance between SNPs in different loci (default: 250000)

Finally, specifiy the number of traits (whose GWAS summary statistics are present at the input path) for which the script should be run in the `#SBATCH --array` option. Array job 


2) Create pan-loci table
To identify pan-loci, that is mega loci including overlapping trait-specific loci, run the **`cntl/p10_locus_lister.sbatch`** script, providing:\

3) Run coloc


## To do
- Patch `coloc.abf()` function to account for sample overlap
- ~~Include in `locus.breaker` function the possibility of providing LOG10 p-value~~ DONE
- Slim down output files? --> Implement new p-value filter discussed at last meeting
- What about COJO collinearity? --> Use tryCatch()
- Include possibility to perform coloc on a specified custom region (skipping locus.breaker and locus.lister steps)
- Include TileDB format as input --> check https://dirk.eddelbuettel.com/papers/useR2021_tiledb_tutorial.pdf
- Include example data to test the pipeline?
- Create a master .sbatch?



## Authors and acknowledgment
Nicola Pirastu ([nicola.pirastu@fht.org](nicola.pirastu@fht.org))\
Arianna Landini ([arianna.landini@external.fht.org](arianna.landini@external.fht.org))\
Sodbo Sharapov ([sodbo.sharapov@fht.org](sodbo.sharapov@fht.org))