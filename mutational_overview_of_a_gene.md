Hello everyone, here is the example code for how to have an overview of alll the mutationsof a particular gene across all the different cancer types in TCGA:

- First we need to download and save the MAF files in your local machine to avoid downloading them everytime.

```R
# Load TCGAbiolinks to access and download TCGA data
library(TCGAbiolinks)
# Load maftools that has functions to manipulate maf files
library(maftools)

# Get all "TCGA" projects stored in GDC
listofprojects <- TCGAbiolinks:::getGDCprojects()$project_id # This gets all the projects in GDC
tcga_projects <- listofprojects[startsWith(listofprojects, "TCGA")] # This selects only project names that start with TCGA
tcga_projects <- sort(tcga_projects) # This sorts them alphabetically

# Create an empty list to store the MAF tibbles
maf_list <- list()
for(i in 1:length(tcga_projects) )
{
  # Read Download the MAF tibble
  maf_tibble <- TCGAbiolinks::GDCquery_Maf(gsub("TCGA-","",tcga_projects[i]), pipelines = "mutect2")
  # Parse it a bit
  colnames(maf_tibble)[colnames(maf_tibble) %in% 'HGVSp_Short'] <- "Protein_Change"
  maf_tibble[,"Tumor_Sample_Barcode"] <- lapply(maf_tibble[,"Tumor_Sample_Barcode"], function(x) substr(x, 1, 12) )
  # Read into maftools
  maf <- read.maf(maf = maf_tibble)
  # Save it into the maf list
  maf_list[[i]] <- maf
}

# Save the maflist into your local machine
saveRDS(maf_list, "~/Desktop/maflist.RDS")

# read it everytime so you avoid having to re-download the data each time
maf_list <- readRDS("~/Desktop/maflist.RDS")

# The results in the list correspond to the names of the projects in the vector "tcga_projects"
# for example the first maf in maf_list corresponds to TCGA-ACC
print(tcga_projects[1])
# "TCGA-ACC"

# You can access each MAF in the list with maf_list[[n]], where n is a number in the list, for example:
oncoplot(maf = maf_list[[1]], top = 10)
```

![Rplot](https://user-images.githubusercontent.com/1195488/133899519-7a36833d-7a0a-4f9f-936a-d0e5302b57b5.png)

- With this file on your computer we can have a look on how are the mutations distributed for a given gene across all cancer types: 

```R
# import the oncomatrix function from the github repository
source("https://raw.githubusercontent.com/PoisonAlien/maftools/master/R/oncomatrix.R")

# read all the 33 maf
maf_list <- readRDS("/astar/r_proj/maflist.RDS")

# Check the mutation frequency of a given gene across all the cohorts:
gene_of_interest <- "TP53"

oncomatrix_table_df <- data.frame()
for(i in 1:length(maf_list)){
  # This line extracts the oncomatrix in character format directly for the gene_of_interest
  oncomatrix <- createOncoMatrix( maf_list[[i]],  gene_of_interest)$oncoMatrix
  if( !is.null(oncomatrix) ){
    oncomatrix_table <- cbind(tcga_projects[i], as.character(oncomatrix) )
    oncomatrix_table_df <- rbind(oncomatrix_table_df, oncomatrix_table)
  }
}

# rename the columnbs of the resulting dataframe
colnames(oncomatrix_table_df) <- c("project","alteration_type")

# re-shape the table to add a count of how many of each alteration types are there per project (cancer type)
oncomatrix_table_df <- aggregate(oncomatrix_table_df, by=list(oncomatrix_table_df$project, oncomatrix_table_df$alteration_type), FUN=length)

# re-name the first 3 columns
colnames(oncomatrix_table_df)[1] <- "project"
colnames(oncomatrix_table_df)[2] <- "alteration_type"
colnames(oncomatrix_table_df)[3] <- "freq"

# just keep the first 3 columns
oncomatrix_table_df <- oncomatrix_table_df[,c(1,2,3)]

# load ggplot2
library(ggplot2)

# Stacked barplot with multiple groups
ggplot(data=oncomatrix_table_df, aes(x=reorder(project, -freq,sum), y=freq, fill=alteration_type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Overview of TP53 mutations in cancer" , y = "Alteration frequency", x = "TCGA projects")

```
![TP53_tcga_mut](https://user-images.githubusercontent.com/1195488/133906603-71d9a033-fa23-4a2d-8107-38fb04bbdf03.png)


- To be completely honest, the block on top is kind of an exercise in futility because we already have the counts for mutation_type for each gene in the MAF files, but I wanted to show you all how you can get the counts directly from the oncomatrix. 
- The pre-computed counts of mutation_type per gene can be accessed with `maf_list[[n]]@gene.summary`.
- Additionally, we can get the total number of samples per project by accessing the clinical data data.frame inside the maf object. The dim() function  `dim(maf_list[[n]]@clinical.data)[1]`
