# Accessing the mutation data

Hi guys, in this short post I will show you how to access the mutation data so you can explore the data with the code that you learnt in class.

### Creating a table (matrix) of gene mutation per samples.

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
#               TCGA-AB-2818        TCGA-AB-2859        TCGA-AB-2895        TCGA-AB-2919        TCGA-AB-2945       
# TP53   "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation"
# PIK3CA "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation" "Missense_Mutation"
# TTN      "Frame_Shift_Ins"   "Frame_Shift_Ins"   "Frame_Shift_Ins"   "Frame_Shift_Ins"   ""  

# The numeric oncomatrix contains the same information; only that it is encoded in numerical form.
# We can extract the numeric oncomatrix with this line:
oncomatrix_num <- oncomatrix$numericMatrix

# Taking a quick look with head() shows the same matrix but with numbers
head(oncomatrix_num)

#          TCGA-AB-2818 TCGA-AB-2859 TCGA-AB-2895 TCGA-AB-2919 TCGA-AB-2945 TCGA-AB-2861 TCGA-AB-2925 TCGA-AB-2931 TCGA-AB-2802
# TP53            1            1            1            1            1            5            1            1            1
# PIK3CA          1            1            1            1            1            0            0            0            0
# TTN             2            2            2            2            0            2            2            2            0

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

With the data in this matrix format, we can start to do some exploration, with some of the functions we saw in class and do more customized analysis:
```R
# We can see the distribution of TP53 mutations in the breast cancer cohort; 
# notice that the first 358 with no label on top means that there are 358 patients with no TP53 mutation.

table( oncomatrix_char["TP53",] )

#                    Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del Missense_Mutation         Multi_Hit 
#              358                46                15                 5               198                 6 
# Nonsense_Mutation       Splice_Site 
#               43                25 

# If we store the results of table() in a variable we can add a name to that column:
# Save the results in a table
tp53_table = table( oncomatrix_char["TP53",] )

# the function names() let's us access or modify the names
names(tp53_table)
# [1] "In_Frame_Del"      "Multi_Hit"         "Frame_Shift_Ins"   "Splice_Site"       "Nonsense_Mutation" "Frame_Shift_Del"  
# [7] "Missense_Mutation" ""   

# Position 8 "" is the one lacking a name; we can set a name like this:
names(tp53_table)[8] <- "No_Mutation"

# We can sort the values by frequency with the function sort()
tp53_table <- sort(tp53_table)

# Lastly if we conver this into data.frame, we can plot it using ggplots as we saw in class.
tp53_df <- as.data.frame(tp53_table)

# Change the column names of the data frame
colnames(tp53_df) <- c("type","freq")

# Plot with ggplots2
library(ggplot2)
p<-ggplot(data=tp53_df, aes(x=type, y=freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() + ggtitle("TP53 mutation distribution")
p
```

![tp53_mutation_dist_brca](https://user-images.githubusercontent.com/1195488/131877413-224af449-e753-4758-ae57-33704dae1a9a.png)

### Other interesting in-built plotting functions.
There are some other plots/analysis you can do very easily because they are already implemented:

##### Ration of transition and transversion:
```R
# Transition and Transversions.
maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = maf.titv)
```
![titv](https://user-images.githubusercontent.com/1195488/131879240-0bcbdd2e-43c9-4ca6-afab-2c398c63dd61.png)


##### Lollipop plot of TP53 mutations in breast cancer tumor:
```R

#lollipop plot for TP53, which is the most frequent mutated gene in the TCGA breast cancer tumors.
lollipopPlot(
  maf = maf,
  gene = 'TP53',
  AACol = 'Protein_Change',
  showMutationRate = TRUE
)
```
![lolipop_tp53](https://user-images.githubusercontent.com/1195488/131879224-1e6fee3f-2a8b-40a7-b1c7-6306847203d8.png)


##### Mutual-exclusity / co-occurrence in mutations
```R
# This plot shows if a particular gene is found to be mutated together with another gene across many tumors (or if their mutations are mutually exclusive)
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
```
![coOcurrance](https://user-images.githubusercontent.com/1195488/131879632-194beca9-5695-4ff0-bfb3-30c49c6c9b6b.png)


##### Druggable targets
```R
# I noticed that sometimes Rstudio plots the figures too small or with weird margins, If this happens reset the the dimensions of the plot with
par(mfrow=c(1,1))

# Plot the druggable targets depending on commonly mutated genes in the tumors of the cohort
# The plot shows biological functions and the druggable targets:
dgi = drugInteractions(maf = maf, fontSize = 1)
```
![druggable](https://user-images.githubusercontent.com/1195488/131880601-d39006a7-1283-4f47-b7e9-14bfc4e0c447.png)

If you want to explore more in detail the information of the available drugs for a given gene, you can format the results like this:
```R
CDH1.dgi = drugInteractions(genes = "CDH1", drugs = TRUE)
CDH1.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

#    Gene interaction_types    drug_name     drug_claim_name
# 1: CDH1                     VOLASERTIB          Volasertib
# 2: CDH1                      ERLOTINIB           Erlotinib
# 3: CDH1                    SELUMETINIB         Selumetinib
# 4: CDH1                   BICALUTAMIDE        Bicalutamide
# 5: CDH1                                                N/A
# 6: CDH1                                PROTEOLYTIC ENZYMES
# 7: CDH1                   CAPECITABINE        Capecitabine
# 8: CDH1                        BI-2536              BI2536
# 9: CDH1                                          ANTISERUM
#10: CDH1                      LAPATINIB           Lapatinib
```
