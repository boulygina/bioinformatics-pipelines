These are some useful instruments for other _in silico_ lab's tasks.

### _E. coli_ strains identification via [StrainEst](https://github.com/compmetagen/strainest)
```bash
$ shh node6  
$ docker run --rm -v /data2:/data2 -v /data1:/data1 -v /data:/data -it compmetagen/strainest bash
$ cd /data2/bio/strainest
$ bowtie2 --very-fast --no-unal -x E_coli//bowtie/align -1 sample_reads_r1.fastq -2 sample_reads_r2.fastq -S sample.sam
$ samtools view -b sample.sam > sample.bam
$ samtools sort sample.bam -o sample.sorted.bam
$ samtools index sample.sorted.bam
$ strainest E_coli/snp_clust.dgrp sample.sorted.bam sample
```

### MLST typing via [srst2](https://github.com/katholt/srst2#basic-usage---mlst)
```bash
$ ssh node6
$ docker run --rm -v /data2:/data2 -v /data1:/data1 -v /data:/data -it estrain/srst2 bash
$ getmlst.py --species "Your_genus your_species"
### Follow the "Suggested srst2 command for use with this MLST database:"
$ srst2 --input_pe sample_reads_r1.fastq sample_reads_r2.fastq \\
--forward _r1 --reverse _r2 --output sample_MLST --log --mlst_db \\
Suggested_genus_species.fasta --mlst_definitions suggested.txt --mlst_delimiter '_'
```
### Evaluating bacterial pathogenicity via [MP3](http://metagenomics.iiserb.ac.in/mp3/)
Run prokka, get the `faa` folder with fasta files, copy them into the `~/scripts/mp3v.1.0' directory. For each fasta file we need to know the minimum amino acid sequence length. For that, run:
```bash
$ cat Bifidobacterium_bifidum_85B.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```
For each genome, run MP3:
```bash
$ ./mp3 file.fasta 1 29 0.2
```
where:
- `1` is for genomic sequence (`2` for metagenome)
- `29` is the minimum amino acid sequence length, which was evaluated previously
- `0.2` is a threshold for SVM module

This tool will predict pathogenic and non-pathogenic proteins based on 3 different algorithms: HMM, SVM and Hybrid.

### Gene abundance calculation via [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
```bash
$ ssh node7
$ docker run --rm -v /data2:/data2 -v /data1:/data1 -v /data:/data --net=host -it mimarkelova/bowtie2:latest
$ bowtie2 -x bowtie_index -1 sample_reads_r1.fastq -2 sample_reads_r2.fastq -S sample.sam
$ /data/subread-1.5.2-Linux-x86_64/bin/featureCounts -p -t CDS -g Name -a ref.gff -o sample.counts.txt sample.sam
```

### vcftools, bcftools  
```bash
$ ssh node9
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/vcftools:latest
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it jweinstk/bcftools_and_tabix:latest
```

### Genome variants annotation via local [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
```
$ ssh node6
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it zhangb1/annovar:latest bash  
$ /table_annovar.pl sample.PASS.vcf /data/reference/annovar/ --buildver hg19 \\
-out sample.PASS.anno -remove -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,exac03,avsnp150,dbnsfp33a,clinvar_20170905,cosmic70 \\
-operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
$ /annotate_variation.pl -out effect_annot -build hg19 \\
annovar/B002.PASS.anno.avinput /data/reference/annovar/
```

### Genome mapping quality assessment via [QualiMap](http://qualimap.bioinfo.cipf.es/)
```bash
$ ssh node6
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it cgwyx/qualimap_conda_docker:latest bash  
$ qualimap bamqc -bam sample.bam -outdir mapQC -outformat pdf --java-mem-size=16G
```

Or via Java GUI:
```bash
$ cd ~/scripts/qualimap_v2.2.1/
$ ./qualimap
```

### Single cell RNA-seq analysis via [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)
```bash
$ node7
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it singlecellportal/cell-ranger-count-2.1.1:latest bash
$ cellranger mkfastq --id single_cell --run 180726_NB501097_0022_AHCHH2BGX5/ --samplesheet SampleSheet.csv 
## (add '--qc' flag for read quality filtering)  
$ cellranger count --id=sample_name --transcriptome=/data/reference/homo_sapiens/10X_Genomics/refdata-cellranger-hg19-1.2.0 --fastqs=/data1/bio/single_cell/outs/fastq_path/single_cell/sample_name --sample=sample_name --expect-cells=2000
```

### SNP linkage disequilibrium calculation
[https://ldlink.nci.nih.gov/](https://ldlink.nci.nih.gov/)

### Run python notebook locally
```bash
$ cd ~/Desktop/notebook
$ ~/anaconda3/bin/jupyter notebook
```
### BLAST 2.2.28+
Local BLAST (nt/nr database):
```bash
$ ssh node0
$ blastn -query file.fasta -db ../../ncbi-blast-2.2.30+/db/nt -out file.out \\
-evalue 0.00001 -outfmt '6 qseqid sseqid pident length evalue sskingdoms stitle' \\ 
-num_threads 4 -num_alignments 1
### Open file.out in table editor, remove duplicates and sort the column with BLAST hits (e.g., column #3):
$ awk '{print $3}' file.out | sort | uniq -c | sort -nr > blast_result.txt
```
Local BLAST with custom database:
```bash
### Building database from file 
$ ssh node0
$ makeblastdb -in file.fasta -dbtype nucl -parse_seqids -out database_name -title "Database_title"
### Running BLASTn for some genes, for example
$ blastn -query genes.fasta -db database_name -out genes.out
### Extracting matched database sequence
$ blastdbcmd -db database_name -entry sequence_name_from_genes_out > genes.found.fa
### Build bed file with database sequence name, start (-1 nt) and end coordinates and extract this exact matching sequence
$ fastaFromBed -fi genes.found.fa -bed coords.bed -fo exact_genes.fasta
```
### Bcl2fastq
```bash
$ ssh node0
$ /usr/local/bin/bcl2fastq -r 20 -d 10 -p 10 -w 8 --no-lane-splitting --min-log-level DEBUG --use-bases-mask Y151,I8,I8,Y151 -o Conversion
```

### One-liners
- Average read length in `.fasta`:

```$ awk '{/>/&&++a||b+=length()}END{print b/a}' sample.fna```

- Read length discribution in `.fastq`:

```$ awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' sample.fastq```

- Convert `.fastq` to `.fasta`:

```$ awk 'NR%4==1{print ">" $0} NR%4==2{print "" $0}' file.fastq > file.fasta```
or 
```$ sed -n '1~4s/^@/>/p;2~4p' file.fastq > file.fasta```

- Split `.fasta` to sequences and their names:

```bash
$ awk 'NR % 2 == 1' file.fasta > names.txt
$ awk 'NR % 2 == 0' file.fasta > sequences.txt
```

- Column number in file:

```$ awk '{print NF}' file | sort -nu | tail -n 1```

- Change java version:

```$ sudo update-alternatives --config java```

- Extract sequence from `.fasta` using the first word in its name:

``` $ perl ~/scripts/extractSequence.pl input.fasta _query_name_ > output.fasta``` 

- Remove `>` in the beginning of a line:

```$ sed 's|^[>]*||' input > output```

- Calculate insert size for paired-end reads:

```$ ~/scripts/bbmap/bbmerge.sh in1=sample_reads_r1.fastq in2=sample_reads_r2.fastq ihist=ihist.txt loose```

- Calculate sample coverage:

```$ ~/scripts/bbmap/mapPacBio.sh in=reads.fastq ref=reference.fasta covstat=output.txt```

- Calculate duplicated lines in file:

```$ sort file.txt | uniq -c | grep 'duplicates_number' ```

