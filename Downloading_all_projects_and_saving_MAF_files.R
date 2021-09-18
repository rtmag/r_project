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
