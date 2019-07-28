# Axiom microarray conversion 
To convert the Axiom genotyping data into the `.plink` format:
1. Edit script `/home/biouser/scripts/GWAS_QC_pipeline/transform_bornmut.php` your input data: modify the number and names of the input files (input example: `/data1/bio/boulygina/ildus/plink/bornmut/genotypes/DF_2017120401_ssp_bestrec_ext_4.txt`), the names of output files, save as `transform_project.php`.
2. Run it using php on node6:
`php /home/biouser/scripts/GWAS_QC_pipeline/transform_bornmut.php`
The output is `project.map` and `project.ped`. 
3. The `map` file lacks chromosome numbers. To add them, edit and run the script `/home/biouser/scripts/GWAS_QC_pipeline/add_chrn.py`. Check if everything is OK in the result `test_project.txt` and rename it as `project.map`.
