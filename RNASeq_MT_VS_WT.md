- Here we are going to check how to do the analysis of mutant vs wt cancer samples.

```R
# set wd
setwd("/astar/r_proj/")

# load libraries
library(TCGAbiolinks)
library("SummarizedExperiment")
library(maftools)
library("DESeq2")
source("https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/oncomatrix.R")

# read UVM RNASEQ
uvm <- readRDS("uvm_rnaseq.rds")

# extract the matrix of RNA-Seq counts (without all the additional clinical, exp, sample data)
uvm_exp <- assay(uvm)

# Parse the Sample names to match the ones in our mutation and clinical data
colnames(uvm_exp) <- substr( colnames(uvm_exp), 1,12)

# Replace the Ensemble id (ENS00...) for gene symbols (TP53, DNMT3A, etc...)
rownames(uvm_exp) <- uvm@rowRanges$external_gene_name

```

- Now what whe want is to subdivide the samples into the one that have the mutation in GNAQ

```R
# Get project names
listofprojects <- TCGAbiolinks:::getGDCprojects()$project_id # This gets all the projects in GDC
tcga_projects <- listofprojects[startsWith(listofprojects, "TCGA")] # This selects only project names that start with TCGA
tcga_projects <- sort(tcga_projects) # This sorts them alphabetically

# load our MAF list
maf_list <- readRDS("/astar/r_proj/maflist.RDS")

# find in which position of the maf_list is UVM
uvm_index <- which(tcga_projects=="TCGA-UVM")

# extract the maf for uvm and remove the bif maf_list
maf_uvm <- maf_list[[uvm_index]]
rm(maf_list)

# get the sample id with mutations in GNAQ
oncomatrix <- createOncoMatrix(maf_uvm, "GNAQ")
mut <- names(which(oncomatrix$numericMatrix["GNAQ",] > 0))

# get sample id without mutations
samples_id <- maf_uvm@clinical.data$Tumor_Sample_Barcode
wt <- samples_id[ !(samples_id %in% mut) ]
```

- Now that we have the samples IDs for patients with and without the mutations in GNAQ, we need to cross this info with the RNA-Seq
```R
# first remove normal samples
# for that we keep in our RNA-Seq matrix only samples that are not normal
uvm_exp_tumors <- uvm_exp[, grep("Normal",uvm$sample_type, invert = TRUE)]

# annotate what samples in the uvm expression tumor matrix are wt and whcih ones are mutant samples 
mutation_annotation <- colnames(uvm_exp_tumors)
mutation_annotation[ mutation_annotation %in% wt] <- "GNAQ_WT"
mutation_annotation[ mutation_annotation %in% mut] <- "GNAQ_MUT"
```
- From here on the analysis is similar to the one observed in WT vs NORMAL samples:

```R
# load deseq2
library("DESeq2")

# create a data dataframe annotating information for the samples
# we can annotate by different fields in the clinical dataset

# get clinical data for UVM
uvm_clinical <- GDCquery_clinic(project = "TCGA-UVM", type = "clinical")

# annotate the samples in uvm_exp_tumors, you can either use the annotations used in the clinical information dataframe we just got in the line above,
# but you need to re-order the rows because they have a different order than in the RNA-SEQ matrix.
# Or you can use the clinical info that was downloaded together with the RNA-Seq matrix and no re-oirder is required .

# here we get the re-order of the clinical info dataframe to match that one in the RNA-seq matrix.
clinicalIndex_reorder <- match(colnames(uvm_exp_tumors),uvm_clinical$submitter_id)

# create the rna-seq matrix annotation dataframe
sampleData <- data.frame( sample = colnames(uvm_exp_tumors), 
                          GNAQ_status = mutation_annotation,
                          gender = uvm_clinical$gender[clinicalIndex_reorder],
                          race = uvm_clinical$race[clinicalIndex_reorder],
                          gender = uvm_clinical$gender[clinicalIndex_reorder],
                          stage = gsub("[ABC]", "", uvm$ajcc_pathologic_stage, perl=TRUE)
                          )


# Create the DDS object
dds <- DESeqDataSetFromMatrix(
  countData = uvm_exp_tumors,
  colData = sampleData,
  design = ~GNAQ_status)

# This command performs the differential expression analysis using the Wald test 
dds_Wald <- DESeq(dds)

# log2 fold changes will represent genes with a higher expression in tumor samples, while negative log2 fold changes
dds_res <- results(dds_Wald, contrast=c("GNAQ_status","GNAQ_MUT","GNAQ_WT") )

# at this step we should also get the matrix with "normalized" expression count, its a matrix just like the one we had before,
# with samples on columns and genes on rows, but the counts are normalized now and can be used for plotting.
dds_vst <- varianceStabilizingTransformation(dds_Wald)

# save table
saveRDS(dds_res,"/astar/r_proj/uvm_gnaq_dds_res.rds")

# save vst matrix
saveRDS(dds_vst,"/astar/r_proj/uvm_gnaq_dds_vst.rds")
```

- We can check the PCA colored by GNAQ mutation status
```R
plotPCA(dds_vst, intgroup=c("GNAQ_status"))
```
![mut_pca](https://user-images.githubusercontent.com/1195488/135437197-8f27d564-923b-40d5-bfed-9f39db88aaf4.png)

- Or check the PCA colored by gender
```R
plotPCA(dds_vst, intgroup=c("gender"))
```
![gender_pca](https://user-images.githubusercontent.com/1195488/135437215-7ab938df-94c1-4ddd-9255-365f6a1b41a6.png)

- Or check the PCA colored by race
```R
plotPCA(dds_vst, intgroup=c("race"))
```
![race_pca](https://user-images.githubusercontent.com/1195488/135437226-9ea22dab-e0ad-4b7c-a151-96ec31a6fc2e.png)

- Or check the PCA colored by stage
```R
plotPCA(dds_vst, intgroup=c("stage"))
```
![stage_pca](https://user-images.githubusercontent.com/1195488/135437238-14e5882e-dc89-4a89-8481-b5f6d4a13e26.png)


- let's check the number of genes that are up and down regulated 

```R
# load ggplot2
library(ggplot2)

# first we plot all the genes in grey color; pch lets you change the type of symbol that gets plotted
plot(dds_res$log2FoldChange,-log10(dds_res$padj),
     xlab="log2 Fold Change",
     ylab="-log10 P-adjusted values",
     col=alpha("grey",.5),
     pch=20)

# we can add some lines to represent the thresholds we used
abline(v=-1,lty = 2,col="grey")
abline(v=1,lty = 2,col="grey")
abline(h=-log10(0.01),lty = 2,col="grey")

# now let's color the genes that pass those thresholds in red
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange) > 1 & dds_res$padj<0.01],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange) > 1 & dds_res$padj<0.01],
       col="red",pch=20)

# lastly, let's add the number of differentially expressed genes to the plot
legend("topright", paste("Mutant:",length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.01))), bty="n") 
legend("topleft", paste("WT:",length(which(dds_res$log2FoldChange<(-1) & dds_res$padj<0.01))), bty="n") 

```

![uvm_mt_wt_volcano](https://user-images.githubusercontent.com/1195488/135441400-17f9354e-4260-4388-95d7-746d1befa03c.png)

```R


```
