# Multiple traits colocalisation and fine mapping

 - Goal of the pipeline
 - Link to paper(s) of interest


## Description
Performs:

1) Loci identification:\
    1.1) Identification of trait-specific loci\
    1.2) Collapsing of across traits overlapping loci into larger genomic regions

2) Association signals​ untangling​\
    2.1) Identification of independent association signals​ by trait\
    2.2) Leave-one-out approach to “clean” single association signal

3) Colocalisation

4) Finemapping of likely causal variant (for each group of colocalasing traits) using coloc posterior probability by SNP


## Installation
Clone this repository


## Quick start
- Give minimum example (example data needed?)

To obtain trait-specific loci, provide the path where your GWAS summary statistics are located as argument to the `--path` option in the `p09_locus_breaker.sbatch`
```
Rscript --vanilla prj_008_multi_coloc_dev/scripts/locus_breaker_wrap.R \
--path "/path/to/GWAS/summary/statistics/directory"
```
Then run specifying the number of traits (29 in this example)
```
sbatch ./prj_008_multi_coloc_dev/cntl/p09_locus_breaker.sbatch --array=1-29
```
To obtain pan loci, aggregating multiple traits overlapping loci in a single locus, provide the path where your trait-specific loci are located as argument to the `--loci_path` option in the `p10_locus_lister.sbatch`

```
Rscript prj_008_multi_coloc_dev/scripts/locus_lister_wrap.R \
--loci_path "/path/to/trait/specific/loci/table/directory"
```
Then run
```
sbatch ./prj_008_multi_coloc_dev/cntl/p10_locus_lister.sbatch
```
Finally, run pairwise colocalisation analysis for all your traits by providing the path where your pan loci table is located as argument to the `--input` option in the `p11_multi_coloc.sbatch`

```
Rscript --vanilla prj_008_multi_coloc_dev/scripts/multi_coloc_wrap.R \
--input "/path/to/overlapping/genomic/regions/table/name_of_your_table.tsv"
```
Then run specifying the number of pan loci (1328 in this example)
```
sbatch ./prj_008_multi_coloc_dev/cntl/p11_multi_coloc.sbatch --array=1-1328%100
```



## Usage

1) Identification of trait-specific loci

To identify the associated loci specific to each trait in a list, run the **`prj_008_multi_coloc_dev/cntl/p09_locus_breaker.sbatch`** script, providing:\
    `--path`: path to directory or multiple (comma separated) directories where GWAS summary statistics for all your traits of interest are located.\
Note this is the only essential option, if not provided the script will throw an error and stop.\
    `--pref`: prefix to the trait name of provided GWAS summary statistics\
    `--suf`: suffix to the trait name of provided GWAS summary statistics\
Prefix and suffix are necessary to extract only the trait name (e.g. "height", "bmi", etc.) from your GWAS summary statistics file names. Of course, this assumes that all files have been named following the same logic.\
    `--out`: path to the (single) directory where output loci tables will be stored. \
    `--chr`: Name of chromosome column in GWAS summary statistics (default: CHROM)\
    `--pos`: Name of genetic position column in GWAS summary statistics (default: GENPOS)\
While all the listed options are quite relevant, they are not strictly essential: for example, if your GWAS summary statistics have been created using the [GAU Regenie pipeline](https://gitlab.fht.org/genome-analysis-unit/nf-pipeline-regenie) and you are happy with having a `"loci_identification"` directory created in the parent directory of your input folder, then you are fine with defaul arguments.\
    `--pvalue`: Name of p-value of effect column in GWAS summary statistics (default: LOG10P)\
This can either be the raw p-value or the log10 transformed one.\ 
    `--sig_pval`: Significant p-value threshold for top hits (default:5e-08)\
    `--limit`: P-value threshold for loci borders (default: 1e-05)\
    `--hole`: Minimum pair-base distance between SNPs in different loci (default: 250000)

Finally, specify the number of traits (whose GWAS summary statistics are present at the input path) for which the script should be run in the `#SBATCH --array` option.



2) Collapsing of across traits overlapping loci into larger genomic regions

To collapse overlapping loci from multiple traits into genomic regions, run the **`prj_008_multi_coloc_dev/cntl/p10_locus_lister.sbatch`** script, providing:\
    `--loci_path`: path to directory where all trait-specific loci table created by the previous step are stored.\
Note this is the only essential option, if not provided the script will throw an error and stop.\
    `--out`: Path to and name of the pan loci file.\



3) Run colocalisation and fine mapping

To perform association signals​ untangling, colocalisation and fine mapping, run the **`prj_008_multi_coloc_dev/cntl/p11_multi_coloc.sbatch`** script, providing:\
    `--input`: path and name of the pan loci table previously created.\
Note this is the only essential option, if not provided the script will throw an error and stop.\
    `--output` "coloc/multi_coloc" \
    `--chr`: Name of the chromosome column in GWAS summary statistics (default: CHROM)\
    `--pos`: Name of the genetic position column in GWAS summary statistics (default: GENPOS)\
    `--rsid`: Name of the rsid column in GWAS summary statistics (default: ID) \
    `--a1`: Name of the effect allele column in GWAS summary statistics (default: ALLELE1)\
    `--a0`: Name of the non effect allele column in GWAS summary statistics (default: "ALLELE0") \
    `--freq`: Name of the effect allele frequency column in GWAS summary statistics (default: "A1FREQ") \
    `--n`: Name of the sample size column in GWAS summary statistics (default: "N") \
    `--effect`: Name of the effect size column in GWAS summary statistics (default: "BETA") \
    `--se`: Name of standard error of effect column in GWAS summary statistics (default: "SE") \
    `--pvalue`: Name of p-value of effect column in GWAS summary statistics (default: "P") \
    `--grch`: Genomic build of GWAS summary statistics (default: 38) \

Finally, specify the number of genomic regions (identified in the previous step) for which the script should be run in the `#SBATCH --array` option.




## To do
- ~~Patch `coloc.abf()` function to account for sample overlap~~ No need according to Sodbo's simulations
- Check munging function to set "essential" info and "optional" info (code will run anyway if only esential info are provided)
- ~~Include in `locus.breaker` function the possibility of providing LOG10 p-value~~ DONE
- ~~Slim down output files? --> Implement new p-value filter discussed at last meeting~~ DONE
- What about COJO collinearity? --> Use tryCatch()
- Include possibility to perform coloc on a specified custom region (skipping locus.breaker and locus.lister steps)
- Include TileDB format as input --> check https://dirk.eddelbuettel.com/papers/useR2021_tiledb_tutorial.pdf --> Linda is taking care of this
- Include example data to test the pipeline?
- Improve naming consistency:
    1) pan.locus and sub_locus --> either use "." or "_"
    2) sub loci are somentimes called "sub_locus" and sometimes "g1"



## Authors and acknowledgment
Nicola Pirastu ([nicola.pirastu@fht.org](nicola.pirastu@fht.org))\
Arianna Landini ([arianna.landini@external.fht.org](arianna.landini@external.fht.org))\
Sodbo Sharapov ([sodbo.sharapov@fht.org](sodbo.sharapov@fht.org))