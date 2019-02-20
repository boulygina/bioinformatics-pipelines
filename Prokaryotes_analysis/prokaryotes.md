# Prokaryotic genome analysis pipeline
The pipeline includes QC of reads, trimming, draft genome and plasmids assembly and annotation, SNP calling and annotation, MLST typing, and orthologs-based phylogenetic tree construction.
The pipeline was built for Illumina paired-end reads. 

### Reads quality check with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
First, we need to perform some quality filtering:
```sh
$ ssh 10.112.112.35
$ ssh node0
$ cd working_directory/reads
$ fastqc sample_R1_001.fastq.gz
$ fastqc sample_R2_001.fastq.gz
```
### Reads trimming with [TRIMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic) v.0.35
```java
$ java -jar ~/scripts/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 \\
sample_R1_001.fastq.gz sample_R2_001.fastq.gz sample_R1_001_trimmed.fastq.gz \\
sample_R1_001_untrimmed.fastq.gz sample_R2_001_trimmed.fastq.gz \\
sample_R2_001_untrimmed.fastq.gz ILLUMINACLIP:adapters.fasta:2:30:10 \\
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
### Reads trimming with [cutadapt](https://cutadapt.readthedocs.io/en/v1.8.1/index.html) v.1.8.1
```bash
$ cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 50 -o sample_R1_001_trimmed.fastq \\
-p sample_R2_001_trimmed.fastq sample_R1_001.fastq.gz sample_R2_001.fastq.gz
```
-  `-a` and `-A` define the sequences which need to be removed from forward and reverse reads according to the FastQC reports
- `-m` defines the minimal read length after trimming.

### Assembly with [SPAdes](http://cab.spbu.ru/software/spades/) v.3.11.1
_Draft genome assembly:_
```
$ ~/scripts/SPAdes-3.10.1-Linux/bin/spades.py --careful -o spades_trimmed \\
-1 sample_R1_001_trimmed.fastq -2 sample_R2_001_trimmed.fastq
```
_Plasmids assembly:_
```
$ ~/scripts/SPAdes-3.10.1-Linux/bin/spades.py --careful -o spades_trimmed_plasmids \\
-1 sample_R1_001_trimmed.fastq -2 sample_R2_001_trimmed.fastq --plasmid
```
SPAdes could be run via docker on node5:
```
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host \\
-it bryce911/spades-3.11.1:latest
$ spades.py ...
```
As a result, in the output `spades_trimmed` directory `contigs.fasta` file is generated.

### Annotation with [prokka](http://www.vicbioinformatics.com/software.prokka.shtml) v.1.12
On the example of _Klebsiella pneumoniae strain Ctrl_:
```bash
$ ~/scripts/prokka-1.12/bin/prokka --cpu 8 --outdir prokka --prefix Kleb \\
--locustag Kleb --genus Klebsiella --species pneumoniae --strain Ctrl \\
spades_trimmed/contigs.fasta
```
As a result, in the output `prokka` directory 12 annotation files are generated. The most important are:
- faa: protein list with amino acid sequences
- ffn: genes list with nucleotide sequences
- gbk: GenBank annotation
- tbl: gene abbreviation list with CDS start and end positions
- tsv: gene abbreviation list with their names, protein products, and EC number.

For rRNA annotation, run on the local machine:
```bash
$ ~/Desktop/soft/barrnap-master/bin/barrnap contigs.fasta
```

### SNP calling with [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/bcftools.html) v.0.19
Mapping reads to the reference genome with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v.2.1.0:
```bash
$ bowtie2-build reference.fasta {prefix}
```
`{prefix}` is a name for the output index files.
```bash
$ bowtie2 -x {prefix} -1 sample_R1_001_trimmed.fastq -2 sample_R2_001_trimmed.fastq -S sample.sam
$ samtools view -bS sample.sam > sample.bam
$ samtools sort sample.bam -o sample_sorted.bam
$ samtools index sample_sorted.bam
$ samtools mpileup -uD -f reference.fasta sample_sorted.bam > sample.bcf
$ bcftools call sample.bcf -mv > sample.vcf
```
### SNP annotation with [SnpEff](http://snpeff.sourceforge.net/) v.3.6
The directory `/data/snpEff/data/Your_organism_genus_and_species/` should contain the annotation file for reference genome in Genbank format `.gbk`. Rename this file to `.gb`. In the end of the configuration file `snpeff.config` add the following line:
`Your_organism_genus_species_strain.genome	:	Your_organism_genus_species_strain`
Run SnpEff:
```java
$ java -jar snpEff.jar build -genbank -v Your_organism_genus_species_strain
$ java -jar snpEff.jar eff -o txt -no-upstream -no-downstream -v Your_organism_genus_species_strain sample.vcf > sample.txt
```
### MLST typing with [srst2](https://github.com/katholt/srst2#basic-usage---mlst)
Run on the local machine:
```
$ getmlst.py --species "Your_genus your_species"

## Follow the "Suggested srst2 command for use with this MLST database:"
$ srst2 --input_pe sample_R1_001_trimmed.fastq sample_R2_001_trimmed.fastq \\
--forward _R1_001 --reverse _R2_001 --output sample_MLST --log --mlst_db \\
Suggested_genus_species.fasta --mlst_definitions suggested.txt --mlst_delimiter '_'
```
### Orthologs-based phylogenetic tree construction with [OrthoMCL DB](http://orthomcl.org/orthomcl/)
Put all the whole-genome fasta files (as one or multiple contigs) in `~/pipes/ortho/fasta` directory. Names of these fasta files should not contain any symbols, except `_`. Remove `~/pipes/ortho/scripts/tmp` directory. Run the analysis:
```bash
$ cd ~/pipes/ortho/scripts
$ ./runProkkaOrtho.sh
$ ./start.sh
```
The most important output files are:
- ortho/ortho_table.txt: a table with the information, in which organisms the certain gene is observed
- results/CoreGenes.fasta: fasta file with whole-genomes alignment. 

Load `CoreGenes.fasta` into FigTree software (on the local machine), choose "Trees" menu and then "PhyML" (maximum likelihood tree construction) with 100-1000 bootstrap iterations. After the tree is built, save it as `.tree` format. Newick-formatted tree could be visualized using the [iTOL](https://itol.embl.de/upload.cgi) or [treeview](http://etetoolkit.org/treeview/) web instruments.

Phylogenetic tree example:
![Phylogenetic_tree](https://github.com/boulygina/bioinformatics-pipelines/blob/master/Prokaryotes_analysis/tree14.png)

