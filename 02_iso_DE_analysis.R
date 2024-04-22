library(limma)

load("/home/projects2/kvs_chunxu_work/Hema_Pro_Iso_DE/hema_iso_gene_mtx.Rdata") ## isofrom_mtx and gene_mtx
meta_data <- read.table("/home/projects2/kvs_chunxu_work/Hema_Pro_Iso_DE/data/sample_annotation.txt", 
                        header = TRUE, stringsAsFactors = FALSE)
sample_ids <- sub(".*_(a[0-9]+)\\..*", "\\1", colnames(isoform_mtx))
population_info <- meta_data$Population[match(sample_ids, meta_data$SampleNr)]
new_colnames <- paste(sample_ids, population_info, sep = "_")
meta_data$SampleID <- paste(meta_data$SampleNr,meta_data$Population,sep = '_')
colnames(isoform_mtx) <- new_colnames
colnames(gene_mtx) <- new_colnames
map_iso_gene <- read_csv(paste0(data_path,'map_iso_gene_vM34.csv'),col_names = c('Isoform','Gene'))

## log2 transformation
isoform_mtx <- log2(isoform_mtx + 1)
gene_mtx <- log2(gene_mtx + 1)

## Differential Gene Abundance
perform_DGA<- function(matrix_data, meta_data, pop1, pop2) {
  samples_to_include <- meta_data$SampleID[meta_data$Population %in% c(pop1, pop2)]
  matrix_subset <- matrix_data[, samples_to_include]
  matrix_subset <- matrix_subset[rowSums(!is.na(matrix_subset)) >= 3, ]
  total_abundance <- colSums(matrix_subset,na.rm = TRUE)
  group <- factor(meta_data$Population[meta_data$SampleID %in% samples_to_include])
  
  design <- model.matrix(~ 0 + group+total_abundance)
  colnames(design) <- gsub("group", "", colnames(design))
  fit <- lmFit(matrix_subset, design)
  
  df.matrix <- makeContrasts(paste(pop1,pop2,sep = '-'), levels = design)
  fit <- contrasts.fit(fit, df.matrix)
  fit <- eBayes(fit)
  
  top_genes <- topTable(fit, number = 20, sort.by = "P")
  return(top_genes)
}
perform_DGA <- function(matrix_data, meta_data, pop1, pop2) {
  # Extract samples from specific groups
  samples_to_include <- meta_data$SampleID[meta_data$Population %in% c(pop1, pop2)]
  matrix_subset <- matrix_data[, samples_to_include]
  
  # Ensure to remove rows containing at least 3 NA
  matrix_subset <- matrix_subset[rowSums(!is.na(matrix_subset)) >= 3, ]
  
  # Calculate total abundance and use it as a covariate
  total_abundance <- colSums(matrix_subset, na.rm = TRUE)
  
  # Create a group factor and ensure that there are only two levels
  group <- factor(meta_data$Population[meta_data$SampleID %in% samples_to_include])
  if (length(levels(group)) != 2) {
    stop("Error: More than two levels found in the group vector. Please specify only two populations.")
  }
  
  # Create a design matrix without intercepts and include total abundance
  design <- model.matrix(~ 0 + group + total_abundance)
  colnames(design) <- gsub("group", "", colnames(design))
  
  # Fit a linear model and use eBay for differential expression analysis
  fit <- lmFit(matrix_subset, design)
  fit <- eBayes(fit)
  
  # Obtain the top 10 genes with significant differences
  top_genes <- topTable(fit, coef = pop2, number = 10, sort.by = "P")
  
  # Return the result of topTable
  return(top_genes)
}

## Perform differential gene usage analysis
perform_DTU <- function(matrix_data, meta_data, pop1, pop2) {
  # Extract samples from specific groups
  samples_to_include <- meta_data$SampleID[meta_data$Population %in% c(pop1, pop2)]
  matrix_subset <- matrix_data[, samples_to_include]
  matrix_subset <- matrix_subset[rowSums(!is.na(matrix_subset)) >= 3, ]
  
  # Calculate total abundance
  total_abundance <- colSums(matrix_subset, na.rm = TRUE)
  
  # Create population factors
  group <- factor(meta_data$Population[meta_data$SampleID %in% samples_to_include])
  if (length(levels(group)) != 2) {
    stop("Error: More than two levels found in the group vector. Please specify only two populations.")
  }
  
  # Create a design matrix without intercepts and include total abundance
  design <- model.matrix(~ 0 + group + total_abundance)
  colnames(design) <- gsub("group", "", colnames(design))
  
  # Fitting linear models
  fit <- lmFit(matrix_subset, design)
  
  # Perform differential splicing analysis
  gene_ids <- map_iso_gene$Gene[match(rownames(matrix_subset), map_iso_gene$Isoform)]
  spliced <- diffSplice(fit,geneid=gene_ids)
  
  # Obtain significantly different splicing events
  top_spliced <- topSplice(spliced, coef = pop2, test = "F", number = 10)
  
  # Return the result of topSplice
  return(top_spliced)
}

## Perform analysis between HSC and MPP
DGA_results <- perform_DGA(gene_mtx,meta_data,pop1 = 'HSC',pop2 = 'MPP_CD150_neg')
DTU_results <- perform_DTU(isoform_mtx,meta_data,pop1 = 'HSC',pop2 = 'MPP_CD150_neg')

######
samples_to_include <- colnames(gene_mtx)
group <- factor(meta_data$Population[meta_data$SampleID %in% samples_to_include])
design <- model.matrix(~ 0 + group + total_abundance)
colnames(design) <- gsub("group", "", colnames(design))
rownames(design)=meta_data$SampleID
total_abundance <- colSums(gene_mtx, na.rm = TRUE)
contrast.matrix<-makeContrasts("HSC-MPP_CD150_neg",levels=design)
##step1
fit <- lmFit(gene_mtx,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  
##step3
tempOutput = topTable(fit2, coef=1, n=10,sort.by = "P")