

```R
library(TCGAbiolinks)

# This library is required to manipulate the format in which the RNA-Seq data is downloaded
# so please install it with
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

# now load it
library("SummarizedExperiment")

################################################################################
# Downloading RNA-Seq data is a little bit differnt than the rest of the datasets
# First we need to set the query parameters such as the TCGA project name,
# data.category, experimental.strategy, and the most important the WORKFLOW.TYPE.
# This workflow.type is the type of data you want to be returned by the GDCdownload
# function. In this case we will download the "raw" expression counts to show you how to do all
# the processing from scratch, but you can also download processed and normalized files.

# Here we are getting the colon adenocarcinoma cancer 

# prepare the query
query_TCGA <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

# download the data
GDCdownload(query = query_TCGA)

# Here we are formatting the data into an S4 class called "RangedSummarizedExperiment". Don't
# Don't worry too much about what this means, the important bit is to know that appart from the expression counts, 
# which is a matrix of samples (columns), genes (rows) and expression counts; we get clinical info, experimental info
# and sample info stored in the slots of the object. This info can be accessed with the "@" operator
coad <- GDCprepare(query_TCGA)

# extract the matrix of RNA-Seq counts (without all the additional clinical, exp, sample data)
coad_exp <- assay(coad)

# Parse the Sample names to match the ones in our mutation and clinical data
colnames(coad_exp) <- substr( colnames(coad_exp), 1,12)

# Replace the Ensemble id (ENS00...) for gene symbols (TP53, DNMT3A, etc...)
rownames(coad_exp) <- coad@rowRanges$external_gene_name

# The interesting thing about this dataset is that we have a mix of normal and tumor samples,
# There are 41 normal samples and 480 tumor samples. We can check that with:
table(coad$sample_type)
#  Metastatic       Primary Tumor     Recurrent Tumor Solid Tissue Normal
#          1                 478                   1                  41

# Let's save that info as it will come in handy in the future, but making only a distintion between 
# normal and tumor samples
sample_type = coad$sample_type

# The grep here finds all the elements in the vector that contain "Normal"; 
# so only elements with the value "Solid Tissue Normal" get picked up.
# For simplicity, let's change that whole name to just "normal"
sample_type[grep("Normal",sample_type)] <- "Normal"

# By passing the option "invert = TRUE", we can get all the elements that don't contain the word "Normal",
# thus selecting all the tumor samples, lets simplify and call them "tumor"
sample_type[grep("Normal", sample_type, invert = TRUE )] <- "Tumor"

# We will be using DESeq2 to perform the differential expression analysis between normal and tumor samples
# so Install the library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# load it
library("DESeq2")

# we need to create a matrix of conditions or covariates we want to contrast when doing the differential Expression Analysis
# In this case let's run the analysis between tumor and normal samples to know what genes are affected in colon cancer
# So, for this case the matrix of covariates will have only the infor about the sample type we saved earlier.
sampleData <- data.frame( sample = colnames(coad_exp), sample_type = sample_type)

# We need to create a DESEQ DATA SET which is an object that combines the "raw" expression counts, the covariate matrix and 
# a "design" formula which is supposed to denote over what fields of the covariate matrix the statistical comparison is meant to be done
# notice the use of "~" symbol
dds <- DESeqDataSetFromMatrix(
  countData = coad_exp,
  colData = sampleData,
  design = ~sample_type)

# This command performs the differential expression analysis using the Wald test 
# means higher expression in normal samples.
dds_Wald <- DESeq(dds)


# With this command we get a table with the results from the differential expression analysis (DEA). Here we are contrasting 
# by sample_type; because we wrote "Tumor" first and "Normal" second in the contrast, positive values of
# log2 fold changes will represent genes with a higher expression in tumor samples, while negative log2 fold changes
dds_res <- results(dds_Wald, contrast=c("sample_type","Tumor","Normal") )

# at this step we should also get the matrix with "normalized" expression count, its a matrix just like the one we had before,
# with samples on columns and genes on rows, but the counts are normalized now and can be used for plotting.
dds_vst <- varianceStabilizingTransformation(dds_Wald)

# Because this might take some time, I recommend you to save the table of results, 
# the normalized count matrix and the original file downloaded from TCGA

# save table
saveRDS(dds_res,"/astar/r_proj/coad_dds_res.rds")

# save vst matrix
saveRDS(dds_vst,"/astar/r_proj/coad_dds_vst.rds")

# save orinal object
saveRDS(coad,"/astar/r_proj/coad_SummarizedExp.rds")

```

- Now that we have the results of the differential expression analysis, let's move to some downstream analysis; including some stuff you saw in class.

```R
# before looking at the results of the DEA, I reccomend to take a look of the PCA
# we can do this easily with the in-bult function to calculate and plot PCA in DESEQ2
# pass the argument of intgroup to color the samples by sample_type
plotPCA(dds_vst, intgroup=c("sample_type"))

# Not surprising the expression profiles of tumor samples are quite different from normal samples.
```

![COAD_RNASEQ_PCA](https://user-images.githubusercontent.com/1195488/134062705-05ec0a8a-dc6d-4d40-85fa-6333ba8454c9.png)

Let's take a look at the table with the DEA results. The table contains:
- the gene name. 
- baseMean; avg. number of read a gene has across all the samples.
- log2FoldChange; the log2 of the ration between the mean number of reads in one conditions / the other conditions; in this case log2( # reads in tumor / # reads in normal)
- lfcSE; the standard error.
- stat; is the Wald statistic (log2FoldChange / lfcSE)
- pvalue; two-tailed pvalue.
- padj; adjusted p-value to account for multiple testing.

```R
head(dds_res)
```
|        |  baseMean |log2FoldChange|lfcSE   |  stat   |  pvalue   |    padj |
|-------|---------|---------------|----------|---------|-----------|---------|
|       | numeric|  numeric   |numeric|numeric|  numeric|numeric|
|TSPAN6  |  5223.831|       0.469351| 0.1448638|  3.239943| 1.19554e-03| 2.77616e-03|
|TNMD    |    42.246|       0.424439| 0.3011284|  1.409495| 1.58689e-01| 2.29012e-01|
|DPM1    |  1747.475|       0.858998| 0.1125037|  7.635284| 2.25322e-14| 1.92205e-13|
|SCYL3   |   534.267|      -0.275480| 0.0695965| -3.958242| 7.55035e-05| 2.16260e-04|
|C1orf112|   359.812|       1.247564| 0.0911147| 13.692232| 1.12984e-42| 6.74517e-41|
|FGR     |   262.740|       0.112187| 0.1940123|  0.578249| 5.63096e-01| 6.49109e-01|


