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

We can get some data. Notice that after playing with the functions for a while, I added some lines to format the data slightly.
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

# The maf tibble is ready to be read by maftools
maf = read.maf(maf = maf_tibble)
```

Now we can use the function we loaded from the internet to generate a matrix of gene mutations (rows) by samples (columns).
```R
# The function createOncoMatrix creates the matrix with rows of the genes we provide (in this case 3 genes) and the samples are in columns. Inside, the content of the matrix indicates if a given sample is mutated in a given gene.
oncomatrix <- createOncoMatrix(maf, c('DNMT3A','FLT3','NPM1'))

# The oncomatrix variable has two matrices with the same information; one written in character and another with numbers encoding the information.
# to get the oncomatrix written in characters you can do this
oncomatrix_char <- oncomatrix$oncoMatrix

# If you check the content with head(); you will see that the matrix denotes the type of mutation a patient has in a given gene.
# Notice that a "" indicates no mutation for that patient in that gene.
head(oncomatrix_char)
#          TCGA-AB-2818        TCGA-AB-2859        TCGA-AB-2895        TCGA-AB-2919        TCGA-AB-2945       
# DNMT3A "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation"
# FLT3   "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation"
# NPM1   "Frame_Shift_Ins"   "Frame_Shift_Ins"   "Frame_Shift_Ins"   "Frame_Shift_Ins"   ""  

# The numeric oncomatrix contains the same information; only that it is encoded in numerical form.
# We can extract the numeric oncomatrix with this line:
oncomatrix_num <- oncomatrix$numericMatrix

# Taking a quick look with head() shows the same matrix but with numbers
head(oncomatrix_num)

#          TCGA-AB-2818 TCGA-AB-2859 TCGA-AB-2895 TCGA-AB-2919 TCGA-AB-2945 TCGA-AB-2861 TCGA-AB-2925 TCGA-AB-2931 TCGA-AB-2802
# DNMT3A            1            1            1            1            1            5            1            1            1
# FLT3              1            1            1            1            1            0            0            0            0
# NPM1              2            2            2            2            0            2            2            2            0

# The code annotation is stored in the oncomatrix variant:
oncomatrix$vc
#                  0                   1                   2                   3                   4                   5 
#                 "" "Missense_Mutation"   "Frame_Shift_Ins"       "Splice_Site"      "In_Frame_Ins" "Nonsense_Mutation" 
#                  6                   7 
#  "Frame_Shift_Del"         "Multi_Hit" 

# 0 = no mutation
# 1 = missesnse mutation
# 2 = Frame shift insertion
# etc ...

```
