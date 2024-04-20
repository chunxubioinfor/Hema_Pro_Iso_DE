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

## log2 transformation and filteration
isoform_mtx <- log2(isoform_mtx + 1)
gene_mtx <- log2(gene_mtx + 1)
#isoform_mtx_filtered <- isoform_mtx[rowSums(!is.na(isoform_mtx)) >= 3, ]
#gene_mtx_filtered <- gene_mtx[rowSums(!is.na(gene_mtx)) >= 3, ]

## Differential Gene Abundance
perform_DE <- function(matrix_data, meta_data, pop1, pop2) {
  samples_to_include <- meta_data$SampleID[meta_data$Population %in% c(pop1, pop2)]
  matrix_subset <- matrix_data[, samples_to_include]
  matrix_subset <- matrix_subset[rowSums(!is.na(matrix_subset)) >= 3, ]
  total_abundance <- colSums(matrix_subset,na.rm = TRUE)
  group <- factor(meta_data$Population[meta_data$SampleID %in% samples_to_include])
  
  if (length(levels(group)) != 2) {
    stop("Error: More than two levels found in the group vector. Please specify only two populations.")
  }
  
  design <- model.matrix(~ 0 + group+ total_abundance)
  colnames(design) <- gsub("group", "", colnames(design))
  
  fit <- lmFit(matrix_subset, design)
  fit <- eBayes(fit)
  
  topTable(fit)
  top_genes <- topTable(fit, coef = pop1, number = 10, sort.by = "P")
}
