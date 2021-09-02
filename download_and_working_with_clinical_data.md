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

# For numerical data we can use summarize; let's check 
# what is the age that patients were diagnosed with cancer

# To accomplish that, we substract the year when they were diagnosed - the year of birth;
# and summarize the results with summary
summary( brca_clinical$year_of_diagnosis - brca_clinical$year_of_birth )
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  26.00   49.00   58.00   58.46   67.00   90.00       4

# It seems like most of the patients were diagnosed with breast cancer around the age of 58.

# another way to do the same exploratory analysis is to take the "age_at_diagnosis" column
# just be aware that the info is in days so we need to divide by 365 to get the years
summary(brca_clinical$age_at_diagnosis/365)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 26.59   49.39   59.07   59.14   68.13   90.06      16 

```

With this clinical data, we can already do lots of exploration like for example if radiation therapy increases the survival of the patients
```R
# As we saw in class today, the column "days_to_last_follow_up" collects the data of 
# when was the last visit of a patient to the clinic after the initial diagnosis.
# we can integrate this data with the information on whether a patient recieved radiation therapy
# to see if we see a diffence between patients that recieved the therapy vs not getting it.

# We can do a plot like we did in class:
library(ggplot2)

ggplot(brca_clinical, aes(y=days_to_last_follow_up, fill=treatments_radiation_treatment_or_therapy)) +
  geom_boxplot() 
```
![radiation](https://user-images.githubusercontent.com/1195488/131886569-74a7d4a9-f57a-48c9-9dd0-9687de3ff42d.png)

The data seems to show that patients that received radiation therapy live slightly longer.

What about drug or chemotherapy?

```R
# The information on drug or chemotherapy is stored in the column: "treatments_pharmaceutical_treatment_or_therapy"
ggplot(brca_clinical, aes(y=days_to_last_follow_up, fill=treatments_pharmaceutical_treatment_or_therapy)) +
  geom_boxplot() 
```
![drugtherapy](https://user-images.githubusercontent.com/1195488/131886698-35df6c3a-5574-44c2-88fa-2e36b973f58f.png)


### Survival analysis 

We can do surival analysis dividing the patients based on gene mutation like we saw in class today. maftools has a function that does this analysis even much easier that the code we saw in class. check it down below:

```R
# Load the maftools library
library(maftools)

# Load the BRCA mutation maf dataset
maf_tibble <- TCGAbiolinks::GDCquery_Maf("BRCA", pipelines = "mutect2")

# Remember to parse the ID of the maf tibble
maf_tibble[,"Tumor_Sample_Barcode"] <- lapply(maf_tibble[,"Tumor_Sample_Barcode"], function(x) substr(x, 1, 12) )

# You might remember that there is a censor variable indicating if we have information on whether the patient
# is alive or deceased; in our clinical data the column is called "vital_status"; with the values "Alive" or "Dead"
# However, the survival analysis function of maftools expects values in binary form:
# where alive=0 and dead=1

# Therefore we need to add a new column into our clinical data frame, were we transform the "Alive" or "Dead" values to 0 or 1.
# To accomplish that we use the code below. I am afraid this is a bit more advanced so please remind me to explain this in one of our upcoming meetings:
vital_status_binary <- brca_clinical$vital_status
vital_status_binary[vital_status_binary == "Alive"] <- 0
vital_status_binary[vital_status_binary == "Dead"] <- 1
vital_status_binary <- as.numeric(vital_status_binary)

# We create a data.frame with all the previous data.frame columns + the new vital_status_binary column
brca_clinical <- data.frame(brca_clinical, vital_status_binary )

# Before we can integrate with the mutation data in the maf file we need to change the name of the column ID to match the same name in the maf file
colnames(brca_clinical)[1] <- "Tumor_Sample_Barcode"

# Now we can integrate the data with the maftools function to read maf
brca_maf = read.maf(maf = maf_tibble, clinicalData = brca_clinical)

# And we can plot survival curves with the following function
par(mfrow=c(1,1)) # you can re-set the ploting parameter with this command if needed 
mafSurvival(maf = brca_maf, 
            genes = 'PIK3CA', 
            time = 'days_to_last_follow_up', 
            Status = 'vital_status_binary', 
)

```
![surv](https://user-images.githubusercontent.com/1195488/131892426-16f9e84a-fd4f-4ba6-a02b-bcb120a9d984.png)

