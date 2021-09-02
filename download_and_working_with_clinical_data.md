# Downloading and working with clinical data

Hi guys, sorry that I was a bit quiet for the past 2 weeks, I was super busy with work but now it is better. 
Today's talk was cool and thought I could show you how to get this kind of data.

```R
# Load TCGAbiolinks to access and download TCGA data
library(TCGAbiolinks)
# Load maftools that has functions to manipulate maf files
library(maftools)

laml_maf_tibble <- TCGAbiolinks::GDCquery_Maf("COAD", pipelines = "mutect2")
colnames(laml_maf_tibble)[colnames(laml_maf_tibble) %in% 'HGVSp_Short'] <- "Protein_Change"

laml_clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
colnames(laml_clinical)[1] <- "Tumor_Sample_Barcode"

laml = read.maf(maf = laml_maf_tibble, clinicalData = laml_clinical)

lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE
)


somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))


mafSurvival(maf = laml, 
            genes = 'TP53', 
            time = 'days_to_last_follow_up', 
            Status = 'vital_status', 
            )



maf_as_dataframe<- head(as.data.frame(laml_maf_tibble))
```