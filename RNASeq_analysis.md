Hi Everyone, in this post I am showing you an example of RNA-Seq analysis for the TCGA dataset of colon adenocarcinoma.

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


Genes with a differential expression between the two conditions (tumor and normal) are usually defined with thresholds on the log2FC and p-adjusted value. In this case let's used the subjective, but stringent threshold of log2fc > 2 (meaning that the gene is expressed 4 times more in tumors compared to normal samples; or 4 times more expressed in normal vs tumor samples), and a strict threshold on padj < 0.001.

```R
# We can find the number or genes with more expression in tumor samples with the following line:
table(dds_res$log2FoldChange>2 & dds_res$padj<0.0001) 
# FALSE  TRUE 
# 51623  3158 

# We can find the number of genes with higher expression in normal samples with this line:
table(dds_res$log2FoldChange<(-2) & dds_res$padj<0.0001) 
# FALSE  TRUE 
# 53304  1477 
```

A popular way to look at the results of a differential expression analysis is with a volcano plot, that combines the log2Fold change and the -log10 P-ajusted values:

```R
# first we plot all the genes in grey color; pch lets you change the type of symbol that gets plotted
plot(dds_res$log2FoldChange,-log10(dds_res$padj),
     xlab="log2 Fold Change",
     ylab="-log10 P-adjusted values",
     col=alpha("grey",.5),
     pch=20)

# we can add some lines to represent the thresholds we used
abline(v=-2,lty = 2,col="grey")
abline(v=2,lty = 2,col="grey")
abline(h=-log10(0.0001),lty = 2,col="grey")

# now let's color the genes that pass those thresholds in red
points(dds_res$log2FoldChange[abs(dds_res$log2FoldChange) > 2 & dds_res$padj<0.0001],
       -log10(dds_res$padj)[abs(dds_res$log2FoldChange) > 2 & dds_res$padj<0.0001],
       col="red",pch=20)
       
# lastly, let's add the number of differentially expressed genes to the plot
legend("topright", paste("Tumor:",length(which(dds_res$log2FoldChange>2 & dds_res$padj<0.0001))), bty="n") 
legend("topleft", paste("Normal:",length(which(dds_res$log2FoldChange<(-2) & dds_res$padj<0.0001))), bty="n") 
```
![coad_volcano](https://user-images.githubusercontent.com/1195488/134067919-5debbcd8-0eec-4390-8266-432a8b4d432e.png)

Another popular way to represent RNA-seq results is with a heatmap:
```R

# We need to install the following libraries

# has interesting computational mathematic methods; 
install.packages("factoextra")

# This library contains the heatmap function Gokce described during class
install.packages("gplots")

# This library allows us to create color palette quickly
install.packages("RColorBrewer")

library(factoextra)
library(RColorBrewer)
library(gplots)

# first we extract from the normalized counts matrix, the genes that are differentially expressed
sig_vst <- dds_vst[ which(abs(dds_res$log2FoldChange) > 2 & dds_res$padj<0.0001) ,]

# If you see, the resulting object is in a "DESeqTransform" format
# we can used the function assay() to extract the values in a matrix format
sig_vst <- assay(sig_vst)

# We can create a color palette, we ask for 9 different colors going from red to blue,
# then we sub-divide those colors into 20 to add granularity and lastly we reverse the order of the colors
# so that highly expressed gener are in red and lowly expressed genes are in blue
colors <- rev(colorRampPalette( brewer.pal(9, "RdBu") )(20))

# Now let's create a color band for the samples based on whether they are tumor or normal:
col_sample_type <- sample_type
col_sample_type[ col_sample_type=="Normal" ] <- "darksalmon"
col_sample_type[ col_sample_type=="Tumor" ] <- "grey5"
  

# This function plots the heatmap of the normalized expression values for the differentially expressed genes.
# notice that the distance function (distfun) which controls the formation of the dendograms uses pearson correlation
# distance measuments from the factoextra library.

# Also, notice that with breaks we are saturating the colors of the heatmap to make the differences more obvious.
# The range of values wiill go from -1.5 to 1.5
# The vector provided to breaks needs to be one element longer than the vector of colors uses.
# In this case 21 because we have 20 colors
# With ColSideColors, we add a color annotation bar on top where normal samples are colored in dark salmon, while tumor samples are in black/grey
heatmap.2(sig_vst, 
          col=colors,
          breaks=seq(-1.5,1.5,length.out=21),
          ColSideColors = col_sample_type,
          scale="row", 
          trace="none",
          distfun = function(x) get_dist(x,method="pearson"),
          labRow = FALSE, labCol = FALSE,
          xlab="Samples", ylab="Genes",
          key.title="Gene expression")
```

![coad_heatmap](https://user-images.githubusercontent.com/1195488/134072567-278296f2-5488-4ff7-a5c9-183b730572b1.png)

Lastly, let's investigate biological pathways that are affected in tumors based on the difference in gene expression.
```R
# Let's install cluster profiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# the previous line should install the packages DOSE and enrichplot directly, if not please install them as shown here
#BiocManager::install("DOSE")
#BiocManager::install("enrichplot")

# lastly get Genome wide annotation database for Human
BiocManager::install("org.Hs.eg.db")


# load all the packages
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("org.Hs.eg.db")


# let's extract the names of genes that are highly expressed in tumors
high_tumor <- rownames(dds_res[ which(dds_res$log2FoldChange>2 & dds_res$padj<0.0001) ,])

# extract the names of genes that are highly expressed in normal samples
high_normal <- rownames(dds_res[ which(dds_res$log2FoldChange<(-2) & dds_res$padj<0.0001) ,])

# this command translates the gene symbol (gene name) of genes with high expression in tumors into it's ENTREZ ID
high_tumor.df <- bitr(high_tumor, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

# do the same for genes with high expression in normal samples
high_normal.df <- bitr(high_normal, fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)

# let's create a list with both, the ENTREZ ID of genes with high expression in tumor and high expression in normal samples.
geneEntrez <- list(high_tumor = high_tumor.df$ENTREZID,
                   high_normal = high_normal.df$ENTREZID
                   )

# run the pathway analysis to indentify the gene onlogies of biological processes (BP) that are affected in tumor samples.
GO_BP <- compareCluster(geneEntrez, fun='enrichGO',
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP")

# Plot the results in a nice dotplot
dotplot(GO_BP, showCategory=15, includeAll=FALSE)
```
![gobp_coad_rnaseq](https://user-images.githubusercontent.com/1195488/134075632-a9e8e891-d686-4bf3-b0f8-30d785207d1a.png)

According to gene onlogies (BP) genes that are highly expressed in tumors are involved in processes of epitelial regulation, this could point out towards an early process of the epithelial-to-mesenchymal transition. Genes that are highly expressed in tumor are involved in immune cell recruitmen, having such genes expressed in cancer cells would make immune system evasion more difficult for tumors.

We 

```R
kegg = compareCluster(geneEntrez, fun="enrichKEGG", organism = "human")

dotplot(kegg, showCategory=15, includeAll=FALSE)
```
![kegg_coad_rna](https://user-images.githubusercontent.com/1195488/134075648-140db84b-dc9f-46ad-bcfe-ce287590ead6.png)

