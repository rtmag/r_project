# Repository to discuss analysis on the R project

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

MAF files are generated after annotating a VCF file. There are many tools to annotate a VCF (vcf2maf, oncotator, funcotator, etc...). Below, I 


|Hugo_Symbol |Entrez_Gene_Id |Center            |NCBI_Build  |Chromosome   |Start_Position  |End_position |Strand  |Variant_Classification  |REF |ALT   |Barcode       |Protein_Change
|-------|------|------------|-------|-------|-------|--------|---------------------------------|--------------|---------------|------------------|-------------|---|
|ABCG4       |64137          |genome.wustl.edu  |37          |11           |119031351       |119031351    |+       |Missense_Mutation SNP   |C   |T     |TCGA-AB-2934  |p.Y567C
|ABL1        |25             |genome.wustl.edu  |37          |9            |133760430       |133760430    |+       |Missense_Mutation SNP   |T   |A     |TCGA-AB-2999  |p.R250W
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417908        |32417909     |+       |Frame_Shift_Ins INS     |-   |ACGG  |TCGA-AB-2839  |p.A170fs       
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417910        |32417911     |+       |Frame_Shift_Ins INS     |-   |ACGG  |TCGA-AB-2844  |p.S169fs
|WT1         |7490           |genome.wustl.edu  |37          |11           |32417909        |32417910     |+       |Frame_Shift_Ins INS     |A   |CGG   |TCGA-AB-2846  |p.S169fs-
|WT1         |7490           |genome.wustl.edu  |37          |11           |32413566        |32413566     |+       |Missense_Mutation SNP   |T   |G     |TCGA-AB-2874  |p.R250W
