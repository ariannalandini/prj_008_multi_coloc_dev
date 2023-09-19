# Multiple traits colocalisation and fine mapping
<br>

## Description
Performs:

1) Loci identification:\
    1.1) Identification of trait-specific loci\
    1.2) Collapsing of across traits overlapping loci into larger genomic regions

2) Association signals​ untangling​\
    2.1) Identification of independent association signals​ by trait\
    2.2) Leave-one-out approach to “clean” single association signal

3) Traits colocalisation

4) Finemapping of likely causal variant (for each group of colocalasing traits) using coloc posterior probability by SNP
<br>


## Installation
Clone this repository
<br>


## Required inputs

**1) GWAS summary statistics**

GWAS summary statistics should use identical column names, irrespective of whether the columns are in a different order or not all columns are present in every set of summary statistics. Additionally, the summary statistics must be based on the same human genome build (either 37 or 38).


| CHROM | GENPOS | ID           | ALLELE0 | ALLELE1 | A1FREQ | N      | BETA    | SE     | LOG10P |
|-------|--------|--------------|---------|---------|--------|--------|---------|--------|--------|
| 1     | 10894	 | 1:10894      | A       | G       | 0.9993 | 367458 |  0.0368 | 0.0763 | 0.2011 |
| 1     | 11171	 | 1:11171      | C       | CCTTG   | 0.9588 | 367458 |  0.0005 | 0.0109 | 0.0161 |
| 1     | 13973	 | 1:13973      | C       | T       | 0.9997 | 367458 | -0.0866 | 0.1608 | 0.229  |
| 1     | 14003	 | 1:14003      | T       | C       | 0.9998 | 367458 | -0.0498 | 0.1627 | 0.1194 |
| 1     | 14022	 | 1:14022      | A       | G       | 0.9998 | 367458 |  0.0461 | 0.1777 | 0.0993 |
| 1     | 23197	 | rs1220638906 | T       | TTAAAA  | 0.9942 | 367458 | -0.038  | 0.0285 | 0.7407 |
<br>

**2) Loci table**

A "guide" table reporting which genomic regions and traits to test for colocalisation. 

the following info (using the below listed column names):

- chr: chromosome number 
- start: starting position of the genomic region of interest
- end: ending position of the genomic region of interest
- trait: name of the associated trait, which will be used in all outputs
- path: path and file name of the GWAS summary statistics file for that trait - **\*.tsv and \*.gz are the file extensions currently accepted**
- type: the type of trait, either "quant" for quantitative traits or "cc" for case-control study\
    if type=="cc"
    - s: for a case-control study, the proportion of samples that are cases (MANDATORY) \
    if type=="quant" and it is known
    - sdY: for a quantitative trait, the population standard deviation of the trait (OPTIONAL - if not given, it will be estimated from the variance of beta and MAF)


| chr | start     | end       | trait    | path                            | type  | sdY  | s    |
|-----|-----------|-----------|----------|---------------------------------|-------|------|------|
| 1   | 11740966  |	11760802  | height   | /group/gwas/height_gwas.cvs.gz  | quant | 0.45 | NA   |
| 1   | 11740966  |	11760802  | bmi      | /random_path/my_gwas_bmi.tsv    | quant | 1.2  | NA   |
| 1   | 17273822  | 17348042  | mpv      | /test/mpv.txt.gz                | quant | 0.02 | NA   |
| 1   | 20560648  | 21242860  | diabetes | /shared/folder/diabetes_103.tsv | cc    | NA   | 0.1  |
| 1   | 28183819  | 28850273  | hbg      | /group/gwas/hgb_gwas.cvs.gz     | quant | 0.01 | NA   |
| 1   | 117583999 | 117640170 | height   | /group/gwas/height_gwas.cvs.gz  | quant | 0.45 | NA   |
| 1   | 117583999 | 117640170 | cvd      | /all_gwas/my_sample_cardio.tsv  | cc    | NA   | 0.36 |
<br>

The loci table can be either produced by using the `p09_locus_breaker.sbatch` and `p10_locus_lister.sbatch` scripts or can be directly provided by the user.

**!!! VERY IMPORTANT !!!**\
When generating the loci table using the `p09_locus_breaker.sbatch` and `p10_locus_lister.sbatch` scripts, it should be noted that the columns labeled as `type`, `sdY`, and/or `s` will not be automatically incorporated into the resulting table, ***unless*** these columns are already present in the GWAS summary statistics file.

In cases where users opt not to make modifications to the GWAS summary statistics, the `type`, `sdY`, and/or `s` columns needs to be manually appended to the table produced by the `p10_locus_lister.sbatch` script. This action should be carried out before running the `p11_multi_coloc.sbatch` script.
<br>

**3) LD reference**
LD reference from UKBB TOPMed imputed genotypes will be used.
<br>


## Quick start

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
<br>


## Options in details

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
<br>


2) Collapsing of across traits overlapping loci into larger genomic regions

To collapse overlapping loci from multiple traits into genomic regions, run the **`prj_008_multi_coloc_dev/cntl/p10_locus_lister.sbatch`** script, providing:\
    `--loci_path`: path to directory where all trait-specific loci table created by the previous step are stored.\
Note this is the only essential option, if not provided the script will throw an error and stop.\
    `--out`: Path to and name of the pan loci file.\
<br>


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
    `--maf`: Minor allele frequency (MAF) thershold removing variants below the set threshold from the analysis (default: 0.0001) \

Finally, specify the number of genomic regions (identified in the previous step) for which the script should be run in the `#SBATCH --array` option.
<br>


## Description of outputs

**Tables:**\
locus_15_final_locus_table.tsv\
Reporting all independent association signals identified for each trait.

locus_15_colocalization.table.all.tsv\
Reporting summary results for all pairwise coloc tests performed.

locus_15_colocalization.table.H4.tsv\
Reporting summary results for pairwise coloc tests performed having PPH4 > 0.9.


**Plots:**\
locus_15_conditioned_loci.pdf\
Visualising the assocation pattern "dissecting" process for each tested trait.

locus_15_pleiotropy_table.pdf\


locus_15_colocalization_plot.pdf\
Regional Manhattan plot of colocalising traits.

locus_15_results_summary_plot.png\
Reporting all independent signals and traits associated at the locus, summarising the colocalising ones and placing them in the context of the chromosome (with genes annotation).
<br>


## References
 - [Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10(5):e1004383](https://doi.org/10.1371/journal.pgen.1004383)
 - [Yang J, Ferreira T, Morris AP et al. . Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat. Genet. 2012;44(4):369–375,S1–S3](doi:10.1038/ng.2213)
<br>


## Authors and acknowledgment
Nicola Pirastu ([nicola.pirastu@fht.org](nicola.pirastu@fht.org))\
Arianna Landini ([arianna.landini@external.fht.org](arianna.landini@external.fht.org))\
Sodbo Sharapov ([sodbo.sharapov@fht.org](sodbo.sharapov@fht.org))
<br>


## To do
- Create conda environment for R packages needed
- Add more possible file extension to the GWAS format in `p09_locus_breaker.sbatch`\
- Possibility to provide costume LD reference (?)\
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