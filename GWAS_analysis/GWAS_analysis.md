# GWAS data analysis
This pipeline describes various types of phenotype association analysis. The PLINK files `clean_project_filter` obtained at the previous QC step are used.

All the proprietary scripts used in this pipeline are located in `~/scripts/GWAS_analysis_pipeline/`.

## Linear association with adjustment for age
Additive model is applied:
```bash
$ ~/scripts/plink_linux_x86_64/plink -bfile clean_project_filter --linear --pheno clean_project_filter.pheno --covar clean_project_filter_age.txt --hide-covar --covar-name age -out output_directory
```
- `clean_project_filter.pheno` is a tab-separated text file with 3 columns: the 1st and the 2nd are sample IDs, the 3rd is a value of a parameter to be analyzed (height, BMI, reaction speed, etc.). For example:

|          |          |     |
| :------- |:---------| :-- |
| sample 1 | sample 1 | 185 |
| sample 2 | sample 2 | 178 |
| sample 3 | sample 3 | 191 |

-  `clean_project_filter_age.txt` is similar to the `pheno` file, but with header and in the 3rd column it contains sample's age:

| FID      | IID      | age |
| :------- |:---------| :-- |
| sample 1 | sample 1 | 14  |
| sample 2 | sample 2 | 12  |
| sample 3 | sample 3 | 12  |

The output is `clean_project_filter_age.assoc.linear`. The format is:

|       |                                                            |
|:----- |:---------------------------------------------------------- | 
|CHR    |Chromosome                                                  |
|SNP    |SNP identifier                                              |
|BP     |Physical position (base-pair)                               |
|A1     |Tested allele (minor allele by default)                     |
|TEST   |Code for the test (see below)                               |
|NMISS  |Number of non-missing individuals included in analysis      |
|BETA/OR|Regression coefficient (--linear) or odds ratio (--logistic)|
|STAT   |Coefficient t-statistic                                     |
|P      |Asymptotic p-value for t-statistic                          |

Since we need the significant association only (p-value ≥ 0.05), discard the other markers:
```bash
$ awk '{if ($NF < 0.05) print $0}' clean_project_filter_age.assoc.linear > clean_project_filter_age.assoc.linear_0.05.txt
```
Add the following headers to the new file:
`CHR	SNP	BP	A1	TEST	NMISS	BETA	STAT	P`

And, finally, assign genes to rsid: save rsid list as `rs.list` and run:
```bash
$ python3 /data2/bio/sandbox/rslist/rslist2flanks.py -i rs.list
```
This script adds the information on the chromosome number, gene name and SNP location to each rsid using the UCSC hg19 database. Output example:

|RS_NAME   |CHROM|POS     |GENE_NAME|PLACE       |
|:-------- |:--- |:------ |:--------|:---------- |
|rs7969362 |chr12|63449964|Y_RNA    |AFTER_100000|
|rs12444268|chr16|20342571|ACSM5    |BEFORE_50000|

## Association with genotype counts results and no adjustment
```bash
$ ~/scripts/plink_linux_x86_64/plink -bfile clean_project_filter --pheno clean_project_filter.pheno -assoc --qt-means -out project_means
$ ~/scripts/plink_linux_x86_64/plink -bfile clean_project_filter --pheno clean_project_filter.pheno -assoc --adjust -out project_adjust
$ ~/scripts/plink_linux_x86_64/plink -bfile clean_project_filter --pheno clean_project_filter.pheno -assoc -out project_beta
```
We have 3 output files as a result:
- `project_means.qassoc.means` with genotype frequencies for each rsid;
- `project_adjust.qassoc.adjusted` with association values;
- `project_beta.qassoc` with beta (regression) coefficients.

From the last file we need to extract SNPs with P-values less than _0.1_ (in this example):
```bash
$ awk '{print $NF,$0}' project_beta.qassoc | sort | cut -f2- -d' ' > project_beta_sorted.qassoc
$ awk '{if ($NF < 0.1) print $0}' project_beta_sorted.qassoc > project_beta_sorted_P_less0.1.qassoc
```
Save these singificant SNPs as single file with rsid:
```bash
$ awk '{print $2}' project_beta_sorted_P_less0.1.qassoc > rs.list
```
Extract genotype frequencies and association values for these SNPs using proprietary scripts (edit them first):
 ```bash
$ nohup ./beta_means.sh > beta_means.log &
$ nohup ./beta_adjust.sh > beta_adjust.log &
$ ~/scripts/plink_linux_x86_64/plink -bfile project_polish_ildus_filter --freq
$ nohup ./freq.sh > freq.log &

$ grep 'GENO' project_beta_means.txt > project_beta_means_GENO.txt
$ grep 'COUNTS' project_beta_means.txt > project_beta_means_COUNTS.txt
$ grep 'FREQ' project_beta_means.txt > project_beta_means_FREQ.txt
$ grep 'MEAN' project_beta_means.txt > project_beta_means_MEAN.txt
$ grep 'SD' project_beta_means.txt > project_beta_means_SD.txt
```
Add the headers: 
- To the `project_beta_sorted_P_less0.1.qassoc`:
` CHR         SNP         BP    NMISS       BETA         SE         R2        T            P `
- To the `project_beta_adjust.txt`:
` CHR         SNP      UNADJ         GC       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY`
- To the `project_freq.txt`:
` CHR         SNP   A1   A2          MAF  NCHROBS`

Save the results in a single tab-separated table in the following order:
- `project_freq.txt`
- `project_beta_sorted_P_less0.1.qassoc`
- `project_beta_adjust.txt`
- `project_beta_means_GENO.txt`
- `project_beta_means_COUNTS.txt`
- `project_beta_means_FREQ.txt`
- `project_beta_means_MEAN.txt`
- `project_beta_means_SD.txt`

Finally, assign genes to rsid: save rsid list as `rs.list` and run:
```bash
$ python3 /data2/bio/sandbox/rslist/rslist2flanks.py -i rs.list
```
Save the result in the 2nd sheet.

