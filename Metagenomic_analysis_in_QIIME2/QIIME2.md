
# Marker gene metagenomics analysis pipeline using QIIME2
The pipeline describes the method of taxonomy identification of bacterial/fungal/eukaryotic community based on the marker gene amplicons (16S/18S rRNA or ITS) on cluster via pre-installed [QIIME](https://docs.qiime2.org/2019.4/) docker, version 2019.4.0.

Note: this particular analysis directory is `/data1/bio/boulygina/qiime2/`.

### Data preprocessing
```bash
ssh node9
docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it qiime2/core
cd working_directory
```
To start the analysis, you should create 2 files: `pe-33-manifest` and `sample-metadata.tsv`.
- The `pe-33-manifest` is a comma-delimited text file which contains reads' absolute filepaths with a header. For example:

sample-id,absolute-filepath,direction
sample1,/data1/bio/boulygina/qiime2/reads/raw/Sample1_R1_001.fastq.gz,forward
sample1,/data1/bio/boulygina/qiime2/reads/raw/Sample1_R2_001.fastq.gz,reverse
sample2,/data1/bio/boulygina/qiime2/reads/raw/Sample2_R1_001.fastq.gz,forward
sample2,/data1/bio/boulygina/qiime2/reads/raw/Sample2_R2_001.fastq.gz,reverse

... etc.

- The  `sample-metadata.tsv` is a tab-delimited text file which contains adapter sequences and samples metadata together with a header and data types. For example:

| #SampleID | BarcodeSequence | LinkerPrimerSequence | Age | Sex |
| :-------- |:----------------| :------------------- | :---- | :---- |
| #q2:types | categorical | categorical | numeric | categorical |
| sample1   | TAAGGCGATAGATCGC |TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG | 35 | m |
| sample2   | GTAGAGGATAGATCGC | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG | 73 | f |
| sample3   | ACTCGCTATAGATCGC | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG | 54 | f |

... etc.

The 16-character "BarcodeSequence" is composed from i7 + i5 Illumina indices used in `SampleSheet.csv`.  "LinkerPrimerSequence" is identical for each sample.
 
### Importing data into QIIME2
```
# import
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path pe-33-manifest --output-path demux.qza --input-format PairedEndFastqManifestPhred33

# statistics
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
```
As a result, 2 files are generated: `demux.qza` which is QIIME2 artifact file (file with data in specific format) and `demux.qzv` which is QIIME2 visualization file. 

`.qzv` files could be viewed via QIIME2  visualization portal: https://view.qiime2.org.
In fact, both `.qza` and `.qzv` files are compressed archives and could be easily  decompressed in terminal or via GUI to extract the necessary information from them.
For more info, refer to https://docs.qiime2.org/2018.8/tutorials/exporting/.

 *=> now we have raw reads number stats (demux.qzv)*

### Sequence QC and feature table construction
Ideally, to set trimming options, you should refer to the `demux.qzv` visualization and set the parameter `--p-trim-left m` (trims off the first m bases of each sequence) and `--p-trunc-len n` (truncates each sequence at position n) according to it. Practically, doing so,  we lose most of the sequences in the text step. This is why I do not set such thresholds and use a bigger number for `--p-max-ee` further while applying DADA2 algorithm: 
```
# trimming, phiX and chimeric reads removal 
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza --p-max-ee 5.0 --p-n-threads 0

# statistics
qiime metadata tabulate --m-input-file stats.qza --o-visualization stats.qzv
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv
```
 *=> now we have filtered reads number stats (stats.qzv), BIOM count table (table.qza), representative sequence variants (rep-seqs.qza) and stats on how many sequences are associated with each sample and with each feature (~OTU) (table.qzv).*
 
 ### Taxonomy profiling
``` 
# mapping of feature IDs to sequences
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

# constructing a phylogenetic tree for further diversity calculations (multiple sequence alignment, masking of highly variable positions, generating a phylogenetic tree, midpoint rooting)
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
```
Taxonomy is assigned to sequences using the Naive Bayes classifier pre-trained on the Greengenes 13_8 99% OTUs; the sequences for `gg-13-8-99-515-806-nb-classifier.qza` have been trimmed to only include 250 bases from the region of the 16S.
Any other classifiers can be found on https://docs.qiime2.org/2019.7/data-resources/.
```
# taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

# statistics
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# generating barplots
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv
```
 *=> now we have bar plots and absolute species abundances (taxa-bar-plots.qzv/data/level-7.csv)*

### Diversity analysis
For exploring species diversity and rarefaction, we need to choose the appropriate sampling depth based on the minimal sequence count in `table.qzv` ("Interactive Sample Detail", on the bottom).
```
# getting all the indices simultaneously
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 4588 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results

# or getting the shannon's index only
qiime diversity alpha --i-table table.qza --p-metric shannon --o-alpha-diversity shannon.qza --output-dir alpha
```
Instead of `shannon` you can use `berger_parker_d`, `simpson_e`, `goods_coverage`,  `pielou_e`, `michaelis_menten_fit`, `chao1`, `lladser_ci`, `chao1_ci`, `esty_ci`, `brillouin_d`, `robbins`, `gini_index`, `fisher_alpha`,  `observed_otus`, `osd`, `menhinick`, `lladser_pe`, `enspie`, `dominance`, `ace`, `doubles`, `singles`, `margalef`, `simpson`, `kempton_taylor_q`, `mcintosh_e`, `strong`, `heip_e`, `mcintosh_d`.

```
# getting the weighted Unifrac metrics
qiime diversity beta-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-metric weighted_unifrac --output-dir beta_weighted

# getting the unweighted Unifrac metrics
qiime diversity beta-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-metric unweighted_unifrac --output-dir beta_unweighted

# getting alpha rarefaction plots
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 10000 --m-metadata-file sample-metadata.tsv --o-visualization alpha_rarefaction.qzv
```
*=> now we have weighted and unweighted Unifracs, alpha rarefaction plots and alpha diversity indices*

For the additional analysis, such as differential abundance testing or PCoA, refer to:
- https://docs.qiime2.org/2019.4/tutorials/moving-pictures/
- https://bioinformaticsworkbook.org/dataAnalysis/Metagenomics/Qiime2.html
