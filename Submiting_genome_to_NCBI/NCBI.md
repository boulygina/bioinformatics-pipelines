# Submitting draft genome to NCBI (in russian)
- Зайти в My NCBI → Submit
- QuickStart → GenBank → Prokaryotic and Eukaryotic Genimes (WGS and complete) → GO
- New Submission

**1. Submitter**

Заполните информацию о вас и о вашей организации.
Звездочкой обозначены обязательные для заполнения поля
→ Continue

**2. General Info**

- Did you already register a BioProject for this research, eg for the submission of the reads to SRA and/or of the genome to GenBank? 
Если ранее не создавали проект для этого генома в NCBI → No

- Did you already register a BioSample for this sample, eg for the submission of the reads to SRA and/or of the genome to GenBank? 
Если ранее не создавали BioSample для этого генома в NCBI → No

- When should this submission be released to the public? 
Выберите «Release immediately following processing», чтобы геном стал доступен сразу после аннотации.  

- Genome info:
-- Assembly method: напр., SPAdes
-- Version or Date program was run: напр., 3.11.1
-- Genome coverage: напр., 790x
-- Sequencing Technology: напр., Illumina MiSeq

- Did your sample include the full genome? → Yes
- Is this the final version? → Yes
- Is it a _de novo_ assembly? → Yes

- Is it an update of existing submission? 
Если геном этого штамма не был ранее загружен в NCBI → No

→ Continue

**3. BioProject General Info**

Public description: краткое описание проекта, для каких целей секвенировали геном.
→ Continue

**4. BioSample Type**

Select the package that best describes your samples: Genome, metagenome or marker sequences (MIxS compliant) → Cultured Bacterial/Archaeal Genomic Sequences MIGS → выбрать подходящую среду обитания

**5. BioSample Attributes**

- Sample Name: название образца
- Organism: напр., Bacillus amyloliquefaciens
- strain: штамм
- collection date: дата сбора образца
- broad-scale environmental context: выбрать biome с сайта, записать ENVO; [пример](https://bioportal.bioontology.org/ontologies/ENVO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01000179)
- local-scale environmental context: выбрать с того же сайта geographic feature; [пример](https://bioportal.bioontology.org/ontologies/ENVO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_00000486)
- environmental medium: выбрать с того же сайта environmental material; [пример](https://bioportal.bioontology.org/ontologies/ENVO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_01001202)
- geographic location: напр., Russia: Kazan
- latitude and longitude: координаты, напр., 38.98 N 77.11 W
- isolation and growth condition: PubMed ID или ссылка для статьи с описанием метода выделения
- number of replicons: 1
- reference for biomaterial: публикация, где впервые упоминается организм

→ Continue

**6. Source**

Annotate this prokaryotic genome in the [NCBI Prokaryotic Annotation Pipeline (PGAP)](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/) before its release → Yes

**7. Files**

- Which of these 3 options describes this genome submission? → 2. One or more chromosomes are still in multiple pieces and/or some sequences are not assembled into chromosomes
- Select file type for the sequences: FASTA
- Select upload type: I will upload all the files now via HTTP/Aspera

Загрузить файл `contigs_ncbi.fasta`
- Do you have AGP files that assemble the individual contigs into scaffolds or chromosomes, OR assemble the submitted gapped sequences into chromosomes? → No

**8. Assignment**
- Is any sequence a complete chromosome? → No
- Does any sequence belong to a plasmid? Если в контигах нет плазмид → No

**9. References**
- Перечислить всех соавторов
- Publication status: Unpublished
- Reference title: примерное название будущей статьи
- Reference authors: Same as sequence authors

**10. Overview**

Проверить, что все в порядке → Submit
