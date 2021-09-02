# Accessing the mutation data

Hi guys, in this short post I will show you how to access the mutation data so you can explore the data with the code that you learnt in class.

##### Creating a table (matrix) of gene mutation per samples.

Let's load the libraries first.
```R
# Load the TCGAbiolinks library
library(TCGAbiolinks)

# Load the maftools library
library(maftools)

# Source the following collection of functions directly from the internet. This allows us to use the function "createOncoMatrix"
source("https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/oncomatrix.R")
```

Now we can get some data. Notice that after playing with the functions for a while, I added some lines to format the data slightly.
```R
# Let's download the breast cancer dataset
maf_tibble <- TCGAbiolinks::GDCquery_Maf("BRCA", pipelines = "mutect2")

# We need to change the name of the column that stores the protein change in the MAF files because some of the maftools functions expect the column
# to be name "Protein_Change", by default the column it is called "HGVSp_Short".
# the data in that column looks like "p.I342V" which means that the mutation is changing the aminoacid 342 in the protein from Isoleucine to Valine.
# the line below changes the name of the column from HGVSp_Short to Protein_Change; There are other ways to rename columns, feel free to explore if you want hehe
colnames(maf_tibble)[colnames(maf_tibble) %in% 'HGVSp_Short'] <- "Protein_Change"

# Another modification is that by default, the data we download comes with the full sample code.
# In TCGA, the IDs look something like this "TCGA-D8-A1XY-01A-11D-A14K-09"; 
# where the "TCGA-D8-A1XY" corresponds to the patient ID, the "01A" whether is tumor or normal tissue,
# and the rest is other codes that indicate if the sample comes from an RNA-Seq, Sequencing or any other experiment.
# In any case, the program expects to have only the patient ID part of the code, thus, the line below removes all the characters after the third "-".
maf_tibble[,"Tumor_Sample_Barcode"] <- lapply(maf_tibble[,"Tumor_Sample_Barcode"], function(x) substr(x, 1, 12) )
```
