

```R
library(survival)
library(survminer)

# Load TCGAbiolinks to access and download TCGA data
library(TCGAbiolinks)
# Load maftools that has functions to manipulate maf files
library(maftools)

# get the oncomatrix function from github
source("https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/oncomatrix.R")

# Get all "TCGA" projects stored in GDC
listofprojects <- TCGAbiolinks:::getGDCprojects()$project_id # This gets all the projects in GDC
tcga_projects <- listofprojects[startsWith(listofprojects, "TCGA")] # This selects only project names that start with TCGA
tcga_projects <- sort(tcga_projects) # This sorts them alphabetically

# read maf list
maf_list <- readRDS("/astar/r_proj/maflist.RDS")

# Threshold on the number of samples that must have the gene mutated to be considered
mutation_thr = .2 # 20% of the samples

surv_pvals_list <- list()

for(i in 1:length(tcga_projects) ){

  ################ Parse survival data  ################ 
  # get clinical data
  clinical <- GDCquery_clinic(project = tcga_projects[i], type = "clinical")
  
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
  
  
  ################ Parse mutation (maf) data  ################ 
  # get sample size
  sample_size <- dim(maf_list[[i]]@clinical.data)[1]
  
  # get genes to consider based on the mutation threshold
  genes2consider <- maf_list[[i]]@gene.summary$Hugo_Symbol[ (maf_list[[i]]@gene.summary$MutatedSamples / sample_size) > mutation_thr ]
  
  # get all the sample ids
  samples_id <- maf_list[[i]]@clinical.data$Tumor_Sample_Barcode
  
  # loop thru all the elements in genes2consider
  pvals_vector <- c()
  for(j in 1:length(genes2consider)){
      oncomatrix <- createOncoMatrix(maf_list[[i]], genes2consider[j])
    
      # samples with  mutations
      mut <- names(which(oncomatrix$numericMatrix[genes2consider[j],] > 0))
    
      # samples with  WT (no mutation)
      wt <- samples_id[ !(samples_id %in% mut) ]
    
      # create a data.frame with the mutation info per patient
      mut_info <- data.frame(sample_id = c(mut , wt ),
                            mut_status = c( rep( "MUT",length(mut)), rep( "WT",length(wt) ) )
      )
    
      ## annotate survival data with whether a patient has a mutation or not in the gene of interest TP53_status
      surv_clinical_gene <- surv_clinical
      surv_clinical_gene$mut_status <- mut_info$mut_status[ match(surv_clinical_gene$bcr_patient_barcode,mut_info$sample_id) ]
    
      # remove NAs from the mut_status
      surv_clinical_gene <- surv_clinical_gene[!is.na(surv_clinical_gene$mut_status),]
    
      ## perform survival analysis
      fit <- survfit(Surv(times, patient.vital_status) ~ mut_status, data = surv_clinical_gene)
      pval <- surv_pvalue(fit)$pval # extract pval
      names(pval) <- genes2consider[j] # add gene name
      pvals_vector <- c( pvals_vector, pval ) # append
  }
  # save pvals_vector
  surv_pvals_list[[i]] <- pvals_vector
  names(surv_pvals_list)[i] <- tcga_projects[i]
}
```

- Check for only significant pvals
```R
# keep only significant pvals
sig_surv_pval <- surv_pvals_list

for(x in 1:length(sig_surv_pval)){
  surv_pvals_list[[x]] <- surv_pvals_list[[x]][surv_pvals_list[[x]]< 0.00001]
  surv_pvals_list[[x]] <- sort(surv_pvals_list[[x]], decreasing = FALSE)
}
```

