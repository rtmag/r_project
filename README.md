# Repository to discuss the analysis for the R project

## Notes on how to access mutation data from TCGA.

Hi everyone. Some of you are familiar with the files that contain the mutation information, but to bring up everyone up to speed let me introduce some file formats. Note, you don't need to know about these file formats to work with them in R but thought it would be interesting for you all to learn. 

#### Variant Call Format (VCF) file:

This is the file that is produced by mutation/variation callers taking as an input the fastq genome sequencing file from WES.
It has one mutation/variant per row and column 10 onwards has the information on whether a given sample in the cohort has that particular mutation.

It looks like this:


|CHROM |POS   |   ID       |  REF   |ALT    |QUAL  |FILTER   |INFO                             |FORMAT       |NA00001         |NA00002          |NA00003
|-------|------|------------|-------|-------|-------|--------|---------------------------------|--------------|---------------|------------------|-------------|
|20     |14370    |rs6054257  |G    | A      |29    |PASS   | NS=3;DP=14;AF=0.5;DB;H2           |GT:GQ:DP:HQ  |0/0:48:1:51,51  |1/0:48:8:51,51   |1/1:43:5:.,.|
|20     |17330    |.          |T     |A     | 3     |q10    | NS=3;DP=11;AF=0.017               |GT:GQ:DP:HQ  |0/0:49:3:58,50  |0/1:3:5:65,3    | 0/0:41:3
|20     |1110696  |rs6040355  |A     |G,T   | 67    |PASS   | NS=2;DP=10;AF=0.333,0.667;AA=T;DB |GT:GQ:DP:HQ  |1/2:21:6:23,27  |2/1:2:0:18,2     |2/2:35:4

- Column 1 tells in which chromosome a mutation is found.
- Column 2 is the position in the genome where the mutation is found.
- If the variant/mutation is known it will be reported in the 3rd column with the ID from the dbSNP database.
- Column 4 shows the nucleotide present in the reference genome.
- Column 5 is the mutation(somatic)/variant(germline) found.
- Column 6 (Qual) is a quality score associated with the algorithm used to find the variants/mutations.
- Column 7 (FILTER) usually tells if there is an issue with the variant/mutation or if there are no issues it says PASS.
- Column 8 (INFO) information regarding the identification of the variant, such as avg number of reads covering the mutation, how prevalent it is in the cohort.
- Column 9 (format) list of fields present in the sample information.
- Column 10 onwards will be the information on whether the mutation/variant is present for each sample in the cohort.

The issue with this format is that we have no information on whether a particular mutation is in a gene, or what gene, or what kind/class of mutation/variation it results in (SNP, insertion, deletion, protein change, etc.). To know this we need to "annotate" the VCF.

#### Mutation Annotation Format (MAF) file:

MAF files are generated after annotating a VCF file. There are many tools to annotate a VCF (vcf2maf, oncotator, funcotator, etc...). Below, I wrote a truncated MAF file with the most relevant columns (it is common to see MAF files with over 100 columns, but you will get the idea of what kind of information is included from this toy example). You can see that it includes information such as gene mutated, type of mutation, protein change and many other that I didn't include here like mutated exon, codon change, etc...

What it is worth noticing is that each row/line corresponds to a mutation/variant found in a particular sample. So if a given gene, let's say WT1 is found mutated in 4 different samples in the cohort, you will have 4 lines with the mutation information (one line per sample).


|Hugo_Symbol |Entrez_Gene_Id |Center            |NCBI_Build  |Chromosome   |Start_Position  |End_position |Strand  |Variant_Classification  |REF |ALT   |Barcode       |Protein_Change
|-------|------|------------|-------|-------|-------|--------|---------------------------------|--------------|---------------|------------------|-------------|---|
|ABCG4       |64137          |genome.wustl.edu  |37          |11           |119031351       |119031351    |+       |Missense_Mutation SNP   |C   |T     |TCGA-AB-2934  |p.Y567C
|ABL1        |25             |genome.wustl.edu  |37          |9            |133760430       |133760430    |+       |Missense_Mutation SNP   |T   |A     |TCGA-AB-2999  |p.R250W
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417908        |32417909     |+       |Frame_Shift_Ins INS     |-   |ACGG  |TCGA-AB-2839  |p.A170fs       
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417910        |32417911     |+       |Frame_Shift_Ins INS     |-   |ACGG  |TCGA-AB-2844  |p.S169fs
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417909        |32417910     |+       |Frame_Shift_Ins INS     |A   |CGG   |TCGA-AB-2846  |p.S169fs-
|WT1         |7490           |genome.wustl.edu  |37          |11           |32413566        |32413566     |+       |Missense_Mutation SNP   |T   |G     |TCGA-AB-2874  |p.R250W


MAF files are the ones we need to obtain to perform many kind of analysis.

## Easy download of TCGA MAF files using R.

TCGA is organized by "cancer type" or "study". This means that we have a cohort for "Breast invasive carcinoma" samples, another for "Ovarian serous cystadenocarcinoma", other for "Prostate adenocarcinoma". All those studies are independent. The full list of studies can be found here together with their abbreviations: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations

Downloading TCGA data is super easy with R, because there are many APIs (which are a set of functions to access the data in a database). Let me show you an example:

```R
# First install the "BiocManager" package
install.packages("BiocManager")

# Load BiocManager
library(BiocManager)

# Now, using BiocManager, install TCGAbiolinks and maftools
BiocManager::install("TCGAbiolinks")
BiocManager::install("maftools")
```

With these two libraries you can already download a MAF file for a particular cancer type and do cool explorations. Let's download the MAF file for TCGA Acute Myeloid Leukemia (abbreviation "LAML") just because is the first one on the list.

```R
# Load the TCGAbiolinks library
library(TCGAbiolinks)

# Using the function "GDCquery_Maf" from the library "TCGAbiolinks", 
# download the MAF file for "LAML" that was created with the WES pipeline "muse" (other possible pipelines are varscan2, somaticsniper and muse).
# !!! NOTE that you have to use the TCGA study abbreviation, in this case "LAML"
maf_tibble <- TCGAbiolinks::GDCquery_Maf("LAML", pipelines = "mutect2")
```

Ee have downloaded and loaded the entire MAF file into a tibble format (don't worry about what tibble is for now, it will be covered in another lesson).

We can now use "maftools" to do some exploration and plots

```R
# Load the maftools library
library(maftools)

# import the MAF tibble into the maftools object
laml = read.maf(maf = maf_tibble)

# generate an oncoprint
oncoplot(maf = laml, top = 10)

# generate other summary plots
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
```
![oncoplot](https://user-images.githubusercontent.com/1195488/129607272-65a79c67-ef0e-4e23-a102-973a3e2f7314.png)
![summaryplot](https://user-images.githubusercontent.com/1195488/129607283-9db3d4dc-3276-49f0-952d-35cd80fbbae1.png)


## Notes on how to check what is inside each "study" in TCGA:

- First you must go to https://portal.gdc.cancer.gov/
- Then on the website type the abbreviation of a cancer study and click on it:
![stp1](https://user-images.githubusercontent.com/1195488/129608691-ab400264-6f8a-4da4-8688-e6d5bf399eb3.PNG)

- In that page you will see how many samples are in that study, how many were sequenced with Whole Exome Sequnecing (WXS) and how many with RNA-Seq:
![stp2](https://user-images.githubusercontent.com/1195488/129608807-50081a2c-48b8-4751-9073-4d7ab35e8533.PNG)


