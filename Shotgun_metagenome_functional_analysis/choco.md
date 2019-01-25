# Shotgun metagenomics analysis pipeline
The pipeline describes the functional annotation of shotgun metagenomic reads using [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home) software and ChocoPhlAn database. 
The pipeline was built for Illumina paired-end reads. 

## Joining reads
The first step is joining the forward and reverse reads together:
```bash
$ cat sample_R1_001.fastq sample_R2_001.fastq > sample1.fastq
```

## Preparing HUMAnN2 environment
We need to replace the demo versions of databases with full versions. For that, check the configuration file after running the HUMAnN2 via docker:
```bash
$ ssh node{7,9}
$ docker run --rm -v /data1:/data1 -v /data2:/data2 -v /data:/data --net=host -it biobakery/humann2:latest bash
$ cd /data2/bio/Chocophlan/
$ humann2_config --print
```
If the config file contains demo versions of databases, download the full versions and update the file:
```bash
$ cd PM
$ humann2_databases --download chocophlan full databases
$ humann2_databases --download uniref uniref90_diamond databases
$ humann2_databases --download utility_mapping full databases
```
If the databases were downloaded once, only the update is needed:
```bash
$ cd PM
$ humann2_config --update database_folders nucleotide databases/chocophlan
$ humann2_config --update database_folders protein databases/uniref
$ humann2_config --update database_folders utility_mapping databases/utility_mapping
```
## Running the mapping
```bash
$ humann2 --input sample1.fastq --output sample1 --threads 30 --gap-fill on --search-mode uniref90 --memory-use maximum
```
The output directory `sample1` contains 3 types of tables:
- `_genefamilies.tsv` contains "gene family abundance reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community. RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples" (- [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)).

(*In Russian*: это число ридов на килобазу, нормализованное на длину гена (RPK). Отражает относительное число копий гена в сообществе. Если образцов несколько, нужно нормализовывать эти значения на глубину покрытия.)
- `_pathabundance.tsv` contains "the abundance of each pathway in the community as a function of the abundances of the pathway's component reactions, with each reaction's abundance computed as the sum over abundances of genes catalyzing the reaction. If the abundances of RXNs 1-4 are [5, 5, 10, 10] in Species_A and [10, 10, 5, 5] in Species_B, HUMAnN2 would report that Species_A and Species_B each contribute 5 complete copies of the pathway. However, at the community level, the reaction totals are [15, 15, 15, 15], and thus HUMAnN2 would report 15 complete copies."

(*In Russian*: это опосредованная (не прямая) сумма чисел всех генов, катализирующих реакции данного метаболизма. Лимитируется числом генов, покрывающих полный путь метаболизма (например, если в пути 1-2-3-4 число генов реакции 1-2 превышает число генов реакции 2-3-4 в 10 раз, то pathabundance будет вычисляться исходя из числа генов реакции 2-3-4.)
- `_pathcoverage.tsv` is "an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance. It is possible for a pathway to be confidently covered at the community level but never confidently detected from any single species".

(*In Russian*: наличие (1) или отсутствие (0) пути метаболизма, неколичественно. 1 = достоверное обнаружение генов всех реакций пути (независимо от его представленности), 0 = менее достоверное обнаружение, т. к. не все гены реакции пути достоверно обнаружены.)

After the mapping, move all the files for each sample to the single directory (e.g. "project").

## Merging and normalization
```bash
$ humann2_join_tables -i project -o project_result/samples_genefamilies.tsv --file_name genefamilies
$ humann2_join_tables -i project -o project_result/samples_pathabundance.tsv --file_name pathabundance
$ humann2_join_tables -i project -o project_result/samples_pathcoverage.tsv --file_name pathcoverage
```
```bash
$ humann2_rename_table --input project_result/samples_genefamilies.tsv --output project_result/samples_genefamilies_names.tsv --names uniref90
$ humann2_renorm_table --input project_result/samples_genefamilies_names.tsv --output project_result/samples_genefamilies_names_cpm.tsv --units cpm --update-snames
$ humann2_renorm_table --input project_result/samples_pathabundance.tsv --output project_result/samples_pathabundance_cpm.tsv --units cpm --update-snames
```
```bash
# Intermediate check, the output should be around 1,000,000.
$ grep -v "#" project_result/samples_genefamilies_names_cpm.tsv | grep -v "|" | cut -f2 | python -c "import sys; print sum(float(l) for l in sys.stdin)"
$ humann2_regroup_table --input project_result/samples_genefamilies_names_cpm.tsv --output project_result/samples_ko_cpm.tsv --groups uniref90_ko
```
Add the following header into the file `samples_pathabundance_cpm.tsv`:
|FEATURE \ SAMPLE | Sample1 | Sample2 | Sample 3 | 
|:----------------|:-------:|:-------:|:--------:|
|     STSite      |  case   |  case   |   ctrl   |
The first header is the samples names, the second - samples group (e.g. "case-control"). Save the file as `samples_pathabund.pcl`.

## Statistical analysis
Explore the association between comparison groups ("case-control" in this case) and samples taxonomy via the Kruskal-Wallis H-test:
```bash
$ humann2_associate --input project_result/samples_pathabund.pcl --last-metadatum STSite --focal-metadatum STSite --focal-type categorical --output project_result/samples_stats.txt
```
The output contains an average pathway abundances in case vs control and the significance of this difference in p- and q-values. 

See which organisms are responsible for the particular pathway: 
```bash
$ grep 'GLYCOGENSYNTH-PWY' samples_pathabund.pcl | less -S
```
Draw a raw barplot for a pathway:
```bash
$ humann2_barplot --input samples_pathabund.pcl --focal-feature  GLYCOGENSYNTH-PWY --focal-metadatum STSite --last-metadatum STSite --output glycogen.png
```
A plot with samples sorted by groups:
```bash
$ humann2_barplot --sort sum metadata --input samples_pathabund.pcl --focal-feature GLYCOGENSYNTH-PWY --focal-metadatum STSite --last-metadatum STSite --output glycogen2.png
```
The output example:
![alt text](https://github.com/boulygina/bioinformatics-pipelines/blob/master/Shotgun_metagenome_functional_analysis/chorismate_biosynthesis_I_sorted.png "Sorted_plot")

Sorting by ecological similarity, normalizing pathway contributions within-sample, and expanding the list of species highlighted:
```bash
$ humann2_barplot --sort similarity --top-strata 12 --scaling normalize --input samples_pathabund.pcl --focal-feature GLYCOGENSYNTH-PWY --focal-metadatum STSite --last-metadatum STSite --output glycogen3.png
```
The output example:
![alt text](https://github.com/boulygina/bioinformatics-pipelines/blob/master/Shotgun_metagenome_functional_analysis/chorismate_biosynthesis_I_grouped.png "Grouped_plot")
Some additional sorting (by the most covered bacteria + by groups + by genera + logarithm):
```bash
humann2_barplot --sort sum metadata -g -a pseudolog --input samples_pathabund.pcl --focal-feature PWY-1269 --focal-metadatum STSite --last-metadatum STSite --output plot.png
```
The association analysis could be done for continuous metadata as well via the spearman correlation (see [HUMAnN2 tutorial](https://bitbucket.org/biobakery/biobakery/wiki/humann2) for details).

