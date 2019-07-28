## Identification of Y haplogroup from NGS data via [Yleaf](https://cluster15.erasmusmc.nl/fmb/Yleaf_v2/index.html?lang=en)
Extract Y chromosome reads from the initial bam:
`samtools view -h sample.mem.sorted.dedup.addrg.bqsr.bam chrY > sample_Y.bam`

On the local machine:
```python
cd /home/lab/Desktop/soft/y_leaf/
python Yleaf.py -bam sample_Y.bam -out out -r 1 -q 20 -b 90 -t 1
```
- `-r` The minimum number of reads for each base above on the quality threshold
- `-q` Minimum quality for each read, integer between 10 and 39, inclusive. If you give it 0, the quality of reads will not be checked
- `-b` The minimum percentage of a base result for acceptance. For example, if you give it 90, then 90% of the reads for each marker should be the same, otherwise that marker will be filtered out
- `-t` Set number of additional threads to use [CPUs] during the alignment process and indexing of BAM files with SAMtools

The output directory is `sample_Y`.

To cite, use `Ralf A., Montiel González D., Zhong K., Kayser M. Yleaf: Software for Human Y-Chromosomal Haplogroup Inference from Next-Generation Sequencing Data. Molecular Biology and Evolution. Volume 35, Issue 5, 1 May 2018, Pages 1291–1294`

## Identification of mtDNA haplogroup from NGS data via [mtDNA-Server](https://mtdna-server.uibk.ac.at/index.html)
Map the reads to the reference mtDNA genome:
```bash
$ ssh node{5-9}
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biocontainers/bwa:v0.7.15_cv2 bash
$ bwa mem -t 24 /data/reference/homo_sapiens/mtDNA/mtDNA_rCRS.fasta sample_R1.fastq.gz sample_R2.fastq.gz > sample_mtDNA.sam
```
Convert sam to bam, sort bam:
```bash
$ ssh node9
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it broadinstitute/gatk:latest bash
$ samtools view -b -S sample_mtDNA.sam > sample_mtDNA.bam
$ samtools sort sample_mtDNA.bam sample_mtDNA.sorted
```
Remove umapped reads:
```
$ samtools view -b -F 4 sample_mtDNA.sorted.bam > sample.sorted.bam
```
Load the result to the mtDNA-Server: https://mtdna-server.uibk.ac.at.
