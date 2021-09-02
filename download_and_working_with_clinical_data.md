# Downloading and working with clinical data

Hi guys, sorry that I was a bit quiet for the past 2 weeks, I was super busy with work but now it is better. 
Today's talk was cool and thought I could show you how to get this kind of data.

```R
# Load TCGAbiolinks to access and download TCGA data
library(TCGAbiolinks)


# To download the clinical data of a particular TCGA project use the function below:
# Notice that here we have to provide the prefix "TCGA-" followed by the name of the project
# Let's use breast cancer in this case
brca_clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# The resulting object is a very large data frame with lots of clinical information for the samples
# you can take a look with head()
head(brca_clinical)

# I am not going to copy the output because it is very large, but you can see that it contains
# info about the stage of the tumor, location, treatment such as radiation or chemo,
# race, age, and the info required for survival analysis.

# As you know we can use the table and summary functions to explore the data a bit

# For categorical data, we can use table() to summarize the dat;
# Let's explore tumor stage, race and gender
table(brca_clinical$ajcc_pathologic_stage)
#   Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III Stage IIIA Stage IIIB Stage IIIC   Stage IV    Stage X 
#        90         86          7          6        358        257          2        155         27         65         20         13

table(brca_clinical$race)
# american indian or alaska native                            asian        black or african american                     not reported 
#                               1                               61                              183                               95 
#                           white 
#                             757 

table(brca_clinical$gender)
# female   male 
#  1085     12 
# Although a bit surprising, you can see men do get breast cancer, only that in low frequency.

# For numerical data we can use summarize; lets check the age of the patients






mafSurvival(maf = laml, 
            genes = 'TP53', 
            time = 'days_to_last_follow_up', 
            Status = 'vital_status', 
            )



maf_as_dataframe<- head(as.data.frame(laml_maf_tibble))
```
