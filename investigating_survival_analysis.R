##############################################################################################################
# SURVIVAL ANALYSIS
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(RTCGA.clinical)

clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
datatable(clinical, filter = 'top', 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
          rownames = FALSE)
##################################################################################################################################
coad_clinical = data.frame(times = clinical$days_to_last_follow_up,
                           bcr_patient_barcode = clinical$bcr_patient_barcode,
                           patient.vital_status = clinical$vital_status)

coad_clinical$patient.vital_status = as.character(coad_clinical$patient.vital_status)

coad_clinical$times[coad_clinical$patient.vital_status == "dead"] <- clinical$days_to_death[coad_clinical$patient.vital_status == "dead"]
# alive=0 and dead=1
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="alive"] = 0
coad_clinical$patient.vital_status[coad_clinical$patient.vital_status=="dead"] = 1

meth.k.id <- data.frame( do.call( rbind, strsplit( names(groups), '-' ) ) )
meth.k.id <- paste(meth.k.id[,1],meth.k.id[,2],meth.k.id[,3],sep="-")

coad_clinical.meth <- coad_clinical[coad_clinical$bcr_patient_barcode %in% meth.k.id,]
coad_clinical.meth <- data.frame(coad_clinical.meth,Methylation.group=groups[match(coad_clinical.meth$bcr_patient_barcode,meth.k.id)])
coad_clinical.meth$Methylation.group <- as.factor(coad_clinical.meth$Methylation.group)
coad_clinical.meth$patient.vital_status <- as.numeric(coad_clinical.meth$patient.vital_status)

#pdf("COAD_methylation_survival.pdf")
kmTCGA(coad_clinical.meth, explanatory.names = "Methylation.group",  pval = TRUE, risk.table=FALSE)
