# Exome analysis pipeline
The pipeline includes single nucleotide polymorphysms (SNP) and copy number variants (CNV) calling and their enhanced annotation (by using freely accessible functional, clinical and gene expression datasets).

All the proprietary scripts used in this pipeline are located in `~/scripts/Exome_analysis_pipeline/`.

## SNP calling and automatic annotation
### Reads quality check with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):
First, we need to perform some quality filtering:
```sh
$ ssh 10.112.112.35
$ ssh node0
$ cd working_directory/reads
$ fastqc sample_R1_001.fastq.gz
$ fastqc sample_R2_001.fastq.gz
```
Before starting to work with files on nodes, modify their permissions with `chmod` (e.g. `sudo chmod 777 working_directory` and `sudo chmod 777 working_directory/*`).

### Reads mapping with [BWA](https://sourceforge.net/projects/bio-bwa/) v.0.7.15 via docker:
For mapping the reads and for all further steps we use hg19 as a human genome reference.
 ```sh
$ ssh node9
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/bwa:latest bash
```
or
```sh
$ ssh node{5-8}
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/bwa:v0.7.15_cv2 bash
```
and then:
```sh
$ cd working_directory
$ bwa mem -t 24 /data/reference/homo_sapiens/ucsc/hg19/sequence/WholeGenomeFasta/genome.fa reads/sample_R1_001.fastq.gz reads/sample_R1_001.fastq.gz > sample.mem.sam
```

### Sorting and indexing with [GATK4](https://software.broadinstitute.org/gatk/gatk4) via docker:
```sh
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it broadinstitute/gatk:latest bash
$ samtools view -b -S sample.mem.sam > sample.mem.bam
$ samtools sort sample.mem.bam sample.mem.sorted
$ samtools flagstat sample.mem.sorted.bam
```
### SNP calling with GATK4:
```javascript
$ gatk --java-options "-Xmx24G" MarkDuplicates -I sample.mem.sorted.bam -O sample.mem.sorted.dedup.bam -M sample_metrics.txt
$ gatk --java-options "-Xmx24G" AddOrReplaceReadGroups -I sample.mem.sorted.dedup.bam -O sample.mem.sorted.dedup.addrg.bam --RGID group1 --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM SureSelect
$ gatk --java-options "-Xmx24G" BaseRecalibrator -R /data/reference/homo_sapiens/ucsc/hg19/sequence/WholeGenomeFasta/genome.fa -I sample.mem.sorted.dedup.addrg.bam -O sample.mem.sorted.dedup.addrg.grp --known-sites /data/reference/homo_sapiens/1000G_omni2.5.hg19.vcf
$ gatk --java-options "-Xmx24G" ApplyBQSR -R /data/reference/homo_sapiens/ucsc/hg19/sequence/WholeGenomeFasta/genome.fa -I sample.mem.sorted.dedup.addrg.bam -bqsr sample.mem.sorted.dedup.addrg.grp -O sample.mem.sorted.dedup.addrg.bqsr.bam
$ gatk --java-options "-Xmx24G" HaplotypeCaller -R /data/reference/homo_sapiens/ucsc/hg19/sequence/WholeGenomeFasta/genome.fa -I sample.mem.sorted.dedup.addrg.bqsr.bam -O sample.GATK.vcf
```
### SNP quality filtering with [vcftools](http://vcftools.sourceforge.net/) v.0.1.x via docker:
 ```sh
$ ssh node9
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/vcftools:latest bash
```
or
```sh
$ ssh node{7,8}
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/vcftools:v0.1.14_cv2 bash
```
and then:
```sh
$ cat sample.GATK.vcf | vcf-annotate -f +/d=20/Q=10/q=10/-a > sample.GATK.filtered.vcf
```
**NOTE**
The bad-quality variants were not filtered yet. Only the column with filter name was added to the file. 

### SNP annotation with ANNOVAR web version:
For ANNOVAR annotation we should compress the file with variants:
```sh
$ gzip sample.GATK.filtered.vcf
```
Then annotate it with [wANNOVAR](http://wannovar.wglab.org/) and download the result named `query.output.exome_summary.txt`.
After that, pick only good variants:
```sh
$ grep 'PASS' query.output.exome_summary.txt > sample_exome_SNP.txt
```

## CNV calling and automatic annotation
### CNV calling with [Canvas](https://github.com/Illumina/canvas) v.1.35:
Rename the file `.bai` to `.bam.bai` in your working directory.
Create the file with heterozygous variants using GATK4: 
```javascript
$ ssh node9
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it broadinstitute/gatk:latest bash
$ gatk VariantFiltration -V sample.GATK.vcf -O sample.GATK_only_het.vcf --genotype-filter-name "not_het" --genotype-filter-expression "(isHomVar == 1)"
$ grep -v 'not_het' sample.GATK_only_het.vcf >  sample.GATK_only_het_filtered.vcf
```
The directory `/data2/bio/CNV_callers/` contains all the necessary files.
Run the CNV calling:
```sh
$ ssh node5
$ cd /data2/bio/CNV_callers/
$ sudo dotnet Canvas-1.35.1.1316.master_x64/Canvas-1.35.1.1316+master_x64/Canvas.dll Germline-WGS -b /working_directory/sample.mem.sorted.dedup.addrg.bqsr.bam -n sample -o sample -r /data/reference/homo_sapiens/ucsc/hg19/sequence/WholeGenomeFasta/genome.fa -g canvas/ref/ -f empty_file.bed --sample-b-allele-vcf /data2/bio/working_directory/sample.GATK_only_het_filtered.vcf
```
The output directory is`/data2/bio/CNV_callers/sample`, and output file is called `CNV.vcf`.

### CNV annotation with [AnnotSV](https://lbgi.fr/AnnotSV/) v.2.0:
Pick only good variants:
```sh
$ grep 'PASS' CNV.vcf > CNV.PASS.vcf
```
Add headers from the initial file `CNV.vcf`.
```sh
$ ssh node0
$ export ANNOTSV="/home/biouser/scripts/AnnotSV_1.1.1"
$ alias AnnotSV='~/scripts/AnnotSV_1.1.1/Sources/AnnotSV-main.tcl'
```
From the working directory:
```sh
$ AnnotSV -SVInputFile CNV.PASS.vcf >& AnnotSV.log &
```
The output directory is `YYYYMMDD_AnnotSV`.

## Annotation of variants using third party databases
Since the automatic annotation we use lacks of full gene nomenclature, their function, expression patterns and other important information requested by collaborators, would be useful to add it using some available databases. We use some proprietary scripts for these purposes.

The directory `/home/biouser/scripts/` contains all the necessary datasets and scripts.

### Gene description
To find gene nomenclature, use [genenames.org](https://www.genenames.org/) service (HUGO Gene Nomenclature Committee).

Copy the list of genes in the annotation obtained earlier and past it into the text editor (e.g. gedit). Replace commas `,` and semicolons `;` with tabs `\t`. Replace the empty rows with "unknown".

Save the first column as `sample_gene_list.txt` and past in to https://www.genenames.org/tools/multi-symbol-checker/. Choose "Approved", "Previous symbols", "Unmatched" and click "Submit". Download the result and merge 2 columns with gene and its nomenclature with previously built database named `gene_description_database.txt`. Replace the empty descriptions (hypothetical gene with names started with "LOC...", etc) with "no". Remove the duplicates and save this new dataset as `gene_description_database.txt`.

Edit and run the script `gene_description_parser.sh`.
The result is `sample_gene_description.txt`.

### Gene function
Gene functions are stored in the mRNA database from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) (table "hgFixed.refSeqSummary"). Genes are presented as transcripts with "NM" ids. So, we need to extract NM_ids from our automatic annotation obtained earlier. For SNP, it could be found in the "AAChange.refgene" column of the annotation file. For CNV, the column named "NM".

Copy and paste this column into the text editor. Replace `\n\n` with `\nUNKNOWN_GENE\n` several times. Save this colums as `sample_NM_list.txt`.

Check `sample_NM_list.txt` for the entries which absent in the database named `summary.tsv` (or in `transcript_ids.txt`).  Add such NMs to the database with 2 "UNKNOWN" columns.

Edit and run the script `gene_function_parser_from_summary.sh`.
The result is `sample_gene_function.txt`.

### Gene expression
Gene expression data were taken from the [GTEx project](https://gtexportal.org/home/) (The Genotype-Tissue Expression). The file with RNA-seq data contains the median TPM values by tissue. 

Check `sample_gene_list.txt` for the entries which absent in the database named `GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct`.  Add such genes to the database with "NA" columns.

Edit and run the script `gene_expression_parser.sh`.
The result is `sample_gene_expression.txt`.

### Gene-drug associations

Information about clinically actionable gene-drug associations and genotype-phenotype relationships is presented in [PharmGKB database](https://www.pharmgkb.org/) (table "Variant and Clinical Annotations Data"). 

Copy the column with RS ids and past it into the table editor. Replace empty rows with dots. Copy the list and paste in into the text editor. Replace dots with "no_rs". Save this list as `sample_rsid_list.txt`.

Check `sample_rsid_list.txt` for the entries which absent in the database named `pharmGKB.txt `.  Add such genes to the database with "unknown" columns.

Edit and run the script `pharmGKB_parser.sh`.
The result is `sample_pharmGKB.txt`.
