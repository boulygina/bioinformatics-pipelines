# GWAS data QC
This pipeline described the quality control procedures for the Genome Wide Association Studies (microarray) data using [PLINK](https://www.cog-genomics.org/plink2) v1.9. The pipeline includes checking for:
- gender mismatch
- missingness for individuals
- missingness for markers
- heterozygosity rate per individual
- relatedness or duplication of individuals
- Hardy-Weinberg equilibrium
- ethnic heterogeneity

As a result, we get the high quality GWAS data in PLINK format.
All the proprietary scripts used in this pipeline are located in `~/scripts/GWAS_QC_pipeline/`.

## Converting microarray data to PLINK format
At first, we need to transform the `.txt` files with whole genome microarray data (e.g. Axiom chip) to PLINK `.ped` and `.map` files. For that, use the script `transform_bornmut.php`. Edit it according your input data: modify the number and names of the input files (input example: `/data1/bio/boulygina/ildus/plink/bornmut/genotypes/DF_2017120401_ssp_bestrec_ext_4.txt`), the names of output files, save as `transform_project.php` and run:
```php
$ ssh node6
$ php transform_project.php
```
The output is `project.ped` and `project.map` files. You need to add chromosome number to `.map`: edit and run `add_chrn.py`. The result is `test_project.map`. Check if everything is OK in this file and rename it as `project.map`.

## Creating and validating binary PLINK files
```bash
$ ssh node0
$ ~/scripts/plink_linux_x86_64/plink --file project --missing-genotype 0 --make-bed --out project

## Output example: 
## Total genotyping rate is 0.988388.
## 786869 variants and 52 people pass filters and QC <-- number of SNPs and individuals in the project.
```
## Applying filters
### At the individual level
***For gender mismatch***

Calculate the mean homozygosity rate across X-chromosome markers for each individual in the study:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile project --check-sex --out project
$ grep 'PROBLEM' project.sexcheck > project.sexprobs
```
Output is `project.sexprobs`. The 3rd column contains gender info in the input file, the 4nd - possible gender according to the array data. Create file with samples with discordant sex information: copy first two columns and save them as *`fail-sexcheck-qc.txt`*. These samples will be filtered out later.

*In Russian:* Для мужчин маркеры на непсевдоаутосомальном участке Х-хромосомы не могут быть гетерозиготными. Такие генотипы обозначаются как "missing". Если много значений "missing" в мужских образцах, то это, скорее всего, женщина. Ожидаемое значение степени гомозиготности для мужчин 1, для женщин < 0.2. Образцы с отклонениями от этих значений удаляются.

***For missing SNPs***

Identification of individuals with elevated missing data rates and outlying heterozygosity rate:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile project --missing --out project
$ ~/scripts/plink_linux_x86_64/plink --bfile project --het --out project
```
In the output `.imiss` the 4th column contains the number of missing SNPs, and the 6th column - the proportion of missing SNPs in individual. The samples with this value between 3-7% need to be filtered out later.

In the output `.het` the 3rd column contains the observed number of heterozygouos genotypes [O(Hom)], and the number of non-missing SNPs in individual [N(NM)]. The heterozygosity rate more or less than the expected value (± 3 SD) may be a sign of contamination and inbreeding, respectively.

To extract such samples, run the R script `imiss-vs-het.r` on the local machine. The plot `project.imiss-vs-het.pdf` will reveal outliers (points outside the red threshold lines). Find their names in output tables sorting the columns "logF_MISS" and "meanHet", and save them as *`fail-imisshet-qc.txt`*.

***For duplicated or related individuals***

To estimate relatedness, we need to calculate the identity by decent, IBD. For each pair of samples the average rate of common alleles is calculated.

For that, markers in high linkage disequilibrium (LD) are removed first to thin the data:
```bash
$ ~/scripts/plink_linux_x86_64/plink --file project --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out project
```
The output file `raw-GWA-data.prune.in` contains SNPs with low LD. Then, generate pairwise IBD for all pairs of individuals in the study based on the reduced marker set:
```bash
$ ~/scripts/plink_linux_x86_64/plink --file project --extract project.prune.in --genome --out project

## Output example: 
## 786869 variants loaded from .bim file.
## 52 people (52 males, 0 females) loaded from .fam.
## --extract: 253831 variants remaining.
## Warning: 4281 het. haploid genotypes present (see test.hh ); many commands treat these as missing.
## Total genotyping rate is 0.992875.
## 253831 variants and 52 people pass filters and QC.
## Excluding 6380 variants on non-autosomes from IBD calculation.
```
`--indep-pairwise 50` means the window size is 50 bp, `0.2` means the SNPs correlation threshold.  IBD is supposed to be 1 for related samples, 0.5 for first-degree relatives, 0.25 for second-degree relatives and 0.125 for third-degree relatives. One sample from each pair with IBD > 0.1875 (halfway between the expected IBD for third- and second-degree relatives) as well as with IBD > 0.98 (considered as duplicates) should be removed. The maximum allowable relatedness in GWAS is second-degree relative.

Identify all pairs of individuals with an IBD > 0.185:
```perl
$ perl run-IBD-QC.pl project
```
Look up the QC-failed samples in *`fail-IBD-QC.txt`*.

***For individuals of divergent ancestry***

Build the PCA with [HapMap3](https://www.sanger.ac.uk/resources/downloads/human/hapmap3.html) samples with known ancestry. "PCA is used to (1) screen the study population for heterogeneous ethnic backgrounds and (2) to correct for potential population stratification (the difference of allele frequencies in ancestral subpopulations)" [from [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007749/)].

The whole procedure is executed as described in the [ENIGMA Network](http://enigma.ini.usc.edu/) [manual](http://enigma.ini.usc.edu/wp-content/uploads/2010/08/ImputationProtocolsv1.0.pdf).

Keep in `project` SNPs from HapMap3 data, save the new dataset as `local`:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile project --missing-genotype 0 --extract HM3.snplist.txt --make-bed --out local
```
Find the strand ambiguous SNPs from the genotyped dataset to avoid strand mismatch among these SNPs:
```bash
$ awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' project.bim | grep -v ambig > local.snplist.txt
```
Remove these SNPs from HapMap3 and save the result as `external`:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile HM3 --extract local.snplist.txt --make-bed --out external
```
Merge `project` and HapMap3 datasets:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile local --bmerge external.bed external.bim external.fam --make-bed --out local_HM3
```
This command will be executed with error: we need to exclude SNPs with multiple alleles from both datasets and then merge them again:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile local --exclude local_HM3-merge.missnp --make-bed --out local_tmp
$ ~/scripts/plink_linux_x86_64/plink --bfile external --exclude local_HM3-merge.missnp --make-bed --out external_tmp
```
Change the values in the last columns in `local_tmp.fam` and `external_tmp.fam`. For `local_tmp.fam`, it should be "3" (if you have only one population in your project; change in to "4", "5" and so on for the next populations). For `external_tmp.fam` it should be ascending numbers (except "9"), starting from "4", for each next population.

Merge new datasets:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile local_tmp --bmerge external_tmp --make-bed --out project_HM3

## Output example: 
## 108142 variants and 1040 people pass filters and QC.
```
Remove `local_tmp.fam` and `external_tmp.fam`.
Plot the PCA in R using `plot-pca-results.r `. Save out the name of outliers in *`fail-PCA-QC.txt`*.

The output example:
![alt text](https://github.com/boulygina/bioinformatics-pipelines/blob/master/GWAS_data_QC/PCA_HapMap3.pdf "PCA_HapMap3")

### At the marker level
Now, combine the samples from all *`fail-QC`* files to one: `project-fail-qc.txt`. Remove these samples from the analysis:
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile project --remove project-fail-qc.txt --make-bed --out clean_project

## Output example: 
## 786869 variants and 48 people pass filters and QC.
```
Apply standard filters for SNPs: Hardy-Weinberg test (0.00001), minor allele frequency value (0.01), genotyping rate (0.05):
```bash
$ ~/scripts/plink_linux_x86_64/plink --bfile clean_project --maf 0.01 --geno 0.05 --hwe 0.00001 --make-bed --out clean_project_filter

## Output example: 
## Total genotyping rate is 0.988795.
## 24751 variants removed due to missing genotype data (--geno).
## --hwe: 1625 variants removed due to Hardy-Weinberg exact test.
## 344008 variants removed due to minor allele threshold(s) (--maf/--max-maf/--mac/--max-mac).
## 416485 variants and 48 people pass filters and QC.
```
The result of this filtering procedure is 5 files in PLINK format (`.bed`, `.bim`, `.fam`, `.hh`, `.log`) ready for further analysis. 

