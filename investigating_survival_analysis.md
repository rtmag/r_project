## Running survival analysis with a different package and with more control over the parsing of the survival information

### Case one GBM PTEN

```R
# Load TCGAbiolinks to access and download TCGA data
library(TCGAbiolinks)
# Load maftools that has functions to manipulate maf files
library(maftools)

# get the oncomatrix function from github
source("https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/oncomatrix.R")

##################################################################################################################
## Get the patients with mutations in a gene ##

# get the MAF tibble
maf_tibble <- TCGAbiolinks::GDCquery_Maf("GBM", pipelines = "mutect2")

# parse the id
maf_tibble[,"Tumor_Sample_Barcode"] <- lapply(maf_tibble[,"Tumor_Sample_Barcode"], function(x) substr(x, 1, 12) )

# read into maftools
gbm_maf = read.maf(maf = maf_tibble)

# oncoplot
oncoplot(maf = gbm_maf, top = 10)

# get all the sample ids
gbm_samples_id <- gbm_maf@clinical.data$Tumor_Sample_Barcode

# get the oncomatrix for PTEN
# Oncomatrix will only return a matrix for patients with mutations 
oncomatrix <- createOncoMatrix(gbm_maf, c("PTEN"))

# samples with PTEN mutations
pten_mut <- names(which(oncomatrix$numericMatrix["PTEN",] > 0))

# samples with PTEN WT (no mutation)
# Oncomatrix will only return a matrix for patients with mutations, so we need to 
# extract, from the total universe of samples, the samples that are not in the pten_mut vector
pten_wt <- gbm_samples_id[ !(gbm_samples_id %in% pten_mut) ]

# create a data.frame with the mutation info per patient
# in this case rep creates a vector of either the word "MUT" or "WT", repeated by the lenght
# of the number of samples with mutation or the # of samples WT
pten_info <- data.frame(sample_id = c(pten_mut , pten_wt ),
                        PTEN_status = c( rep( "MUT",length(pten_mut)), rep( "WT",length(pten_wt) ) )
                        )

##################################################################################################################
## Parse survival data ##

# get clinical data
clinical <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")

# parse survival data
surv_clinical = data.frame(times = clinical$days_to_last_follow_up,
                           bcr_patient_barcode = clinical$bcr_patient_barcode,
                           patient.vital_status = clinical$vital_status)

# re-name the vital status column
surv_clinical$patient.vital_status = as.character(surv_clinical$patient.vital_status)

# The survival times for alive patient is in the column "time", while the times for patients that passed away
# is stored in a different column called "days_to_death"
dead_ix <- which(surv_clinical$patient.vital_status == "Dead") # index to the samples that passed away
surv_clinical$times[dead_ix] <- clinical$days_to_death[dead_ix]

# let's remove NAs and not reported
surv_clinical <- surv_clinical[!is.na(surv_clinical$times),]
surv_clinical <- surv_clinical[surv_clinical$patient.vital_status != "Not Reported",]

# change the vital status to numeric values alive=0 and dead=1
surv_clinical$patient.vital_status[surv_clinical$patient.vital_status=="Alive"] = 0
surv_clinical$patient.vital_status[surv_clinical$patient.vital_status=="Dead"] = 1

# finally let's set patient.vital_status to numeric
surv_clinical$patient.vital_status <- as.numeric(surv_clinical$patient.vital_status)


##################################################################################################################
## annotate survival data with whether a patient has a mutation or not in the gene of interest
surv_clinical$PTEN_status <- pten_info$PTEN_status[ match(surv_clinical$bcr_patient_barcode,pten_info$sample_id) ]

# remove NAs from the PTEN_status
surv_clinical <- surv_clinical[!is.na(surv_clinical$PTEN_status),]

##################################################################################################################
# install RTCGA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RTCGA")

# load RTCGA
library(RTCGA)

## Run survival analysis
kmTCGA(surv_clinical, explanatory.names = "PTEN_status",  pval = TRUE, risk.table=TRUE)

```

![GBM_PTEN](https://user-images.githubusercontent.com/1195488/134184313-cd4c7919-cfbd-4b99-9ef8-61d0e2dd4333.png)



### let's try with TP53

```R
######################################################################################################################################################################################
# get the oncomatrix for TP53
# Oncomatrix will only return a matrix for patients with mutations 
oncomatrix <- createOncoMatrix(gbm_maf, c("TP53"))

# samples with PTEN mutations
tp53_mut <- names(which(oncomatrix$numericMatrix["TP53",] > 0))

# samples with PTEN WT (no mutation)
# Oncomatrix will only return a matrix for patients with mutations, so we need to 
# extract, from the total universe of samples, the samples that are not in the pten_mut vector
tp53_wt <- gbm_samples_id[ !(gbm_samples_id %in% tp53_mut) ]

# create a data.frame with the mutation info per patient
# in this case rep creates a vector of either the word "MUT" or "WT", repeated by the lenght
# of the number of samples with mutation or the # of samples WT
tp53_info <- data.frame(sample_id = c(tp53_mut , tp53_wt ),
                        TP53_status = c( rep( "MUT",length(tp53_mut)), rep( "WT",length(tp53_wt) ) )
)


# parse survival data
surv_clinical = data.frame(times = clinical$days_to_last_follow_up,
                           bcr_patient_barcode = clinical$bcr_patient_barcode,
                           patient.vital_status = clinical$vital_status)

# re-name the vital status column
surv_clinical$patient.vital_status = as.character(surv_clinical$patient.vital_status)

# The survival times for alive patient is in the column "time", while the times for patients that passed away
# is stored in a different column called "days_to_death"
dead_ix <- which(surv_clinical$patient.vital_status == "Dead") # index to the samples that passed away
surv_clinical$times[dead_ix] <- clinical$days_to_death[dead_ix]

# let's remove NAs and not reported
surv_clinical <- surv_clinical[!is.na(surv_clinical$times),]
surv_clinical <- surv_clinical[surv_clinical$patient.vital_status != "Not Reported",]

# change the vital status to numeric values alive=0 and dead=1
surv_clinical$patient.vital_status[surv_clinical$patient.vital_status=="Alive"] = 0
surv_clinical$patient.vital_status[surv_clinical$patient.vital_status=="Dead"] = 1

# finally let's set patient.vital_status to numeric
surv_clinical$patient.vital_status <- as.numeric(surv_clinical$patient.vital_status)


##################################################################################################################
## annotate survival data with whether a patient has a mutation or not in the gene of interest
surv_clinical$TP53_status <- tp53_info$TP53_status[ match(surv_clinical$bcr_patient_barcode,tp53_info$sample_id) ]

# remove NAs from the PTEN_status
surv_clinical <- surv_clinical[!is.na(surv_clinical$TP53_status),]

#########################################################################

## Run survival analysis
kmTCGA(surv_clinical, explanatory.names = "TP53_status",  pval = TRUE, risk.table=TRUE)
```

![p53_gbm](https://user-images.githubusercontent.com/1195488/134184294-5931e0f0-9a5b-4bc2-81c6-c97b4cf13b73.png)
