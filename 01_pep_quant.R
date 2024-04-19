library(parallel)
#Initialization
data_path <- '/home/projects2/kvs_chunxu_work/Hema_Pro_Iso_DE/data/'
work_path <- '/home/projects2/kvs_chunxu_work/Hema_Pro_Iso_DE/'
output_path <- '/home/projects2/kvs_chunxu_work/Hema_Pro_Iso_DE/output/'
mapPeptides <- function(peptides, reference) {
  
  ### Test input
  peptides <- as.character(peptides)
  stopifnot(is.character(reference))
  
  ### Helper function
  mfUnlistAndExtract <- function(aList, pattern) {
    tmp <- unlist(aList)
    tmp <- tmp[str_detect(names(tmp), pattern = pattern)]
    return(tmp)
  }
  
  ### Map and format
  AhoCorasickTrie::AhoCorasickSearch(
    keywords = peptides,
    text = reference,
    alphabet = "aminoacid"
  ) %>% 
    enframe(
      name = 'protein',
      value = 'list'
    ) %>% 
    rowwise() %>% 
    reframe(
      peptide  = mfUnlistAndExtract(list, pattern ='Keyword'),
      protein = protein[1]
      #position = mfUnlistAndExtract(list, pattern ='Offset'),
    ) %>% 
    return()
} ##Helper Function
get_sequence <- function(PrecursorIds){
  PeptideSequence = str_extract(PrecursorIds, "(?<=_).*(?=_.)")
  PeptideSequence_without_mod <- str_remove_all(PeptideSequence, "\\[.*?\\]")
  return(PeptideSequence_without_mod)
}
numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, c("data_path", "output_path", 'work_path')) 
clusterEvalQ(cl, {
  library(tidyverse)
  library(IsoBayes)
  library(AhoCorasickTrie)
  library(parallel)
  library(mailR)
  library(flowTraceR)
})

if(FALSE){
  ### Load annotation data
  gencodeAa <- Biostrings::readAAStringSet(
    paste0(data_path,'gencode.vM34.pc_translations.fa.gz')
  )
  names(gencodeAa) <- map_chr(
    str_split(names(gencodeAa), pattern = '\\|', n = 3),
    ~ .x[2]
  )
  gencodeAa <- unique(gencodeAa)
  
  ### Handle non-ref amino acids
  #gencodeAaTrim <- trimLRPatterns(subject = gencodeAa, Lpattern = 'X', Rpattern = 'X')
  #letterCount <- letterFrequency(gencodeAaTrim, letters = names(AMINO_ACID_CODE))
  #gencodeAaTrim <- gencodeAaTrim[which( rowSums( letterCount[,c('U','O','B','J','Z','X')]) == 0 )]
  
  
  ### Annoation
  gencodeGTF <- rtracklayer::import(paste0(data_path,'gencode.vM34.chr_patch_hapl_scaff.annotation.gtf.gz'))
  gencode44annotation <- as.data.frame(GenomicRanges::mcols(gencodeGTF[which(
    gencodeGTF$type == 'transcript'
  ),]))
  
  gen44pc <- 
    gencode44annotation %>% 
    filter(
      gene_type == 'protein_coding',
      transcript_type == 'protein_coding'
    ) %>% 
    select(gene_id, gene_name, transcript_id, protein_id, transcript_type) %>% 
    distinct()
  
  save(gencodeAa, gencode44annotation, gen44pc, file= paste0(data_path,'mouse_gencode.vM34.protein.Rdata'))
  
}
load(file = paste0(data_path,'mouse_gencode.vM34.protein.Rdata')) # gencodeAa, gencode44annotation, gen44pc
gencodeAaChr <- as.character(gencodeAa)
if(TRUE) {
  library(tidyverse)
  ### Read into the file
  peptide_mtx <- read_tsv(paste0(data_path,'PeptideMatrix.tsv')) %>%
    mutate(PeptideSequence = get_sequence(EG.PrecursorId)) %>%
    select(-PG.ProteinGroups, -PG.GroupLabel,-EG.PrecursorId) %>%
    mutate(across(-PeptideSequence, ~ifelse(. == "Filtered", 0, .)),  
           across(-PeptideSequence, ~as.numeric(gsub(",", ".", .x, fixed = TRUE))))

  sample_names <- colnames(peptide_mtx)[-69]
  PepTable <- 
    peptide_mtx %>%
    dplyr::rename(peptide = PeptideSequence)
  
  ### Map peptides
  PepInfo <- 
    mapPeptides(
      peptides = PepTable$peptide,
      reference = gencodeAaChr
    ) %>% 
    dplyr::rename(
      iso_id = protein
    ) %>% 
    mutate(
      gene_id = gen44pc$gene_id[match(
        iso_id, gen44pc$transcript_id
      )]
    ) %>% 
    group_by(peptide) %>% 
    summarise(
      n_iso  = n_distinct(iso_id),
      n_gene = n_distinct(gene_id),
      iso_ids  = paste(iso_id, collapse = "|"),
      gene_ids = paste(gene_id, collapse = "|"),
      .groups = 'drop'
    ) %>% 
    mutate(
      iso_specific  = n_iso == 1,
      gene_specific = n_gene == 1,
    )
  
  ### Join annotation
  PepMapped <- inner_join(
    PepTable,
    PepInfo,
    by = 'peptide'
  ) %>% 
    select(
      peptide,
      iso_specific, gene_specific, n_iso, n_gene, iso_ids, gene_ids,
      everything()
    )
  
  save(PepMapped, file = paste0(data_path,'mouse_pep_mapped.Rdata'))
}
if(FALSE) {
  load(file = paste0(data_path,'mouse_gencode.vM34.protein.Rdata'))# gencodeAa, gencode44annotation, gen44pc
  head(gen44pc$gene_id)
  bind_rows(gencode44annotation) %>% 
    select(transcript_id, gene_id) %>% 
    distinct() %>% 
    filter(
      transcript_id %in% c(names(gencodeAa))
    ) %>% 
    as_tibble() %>% 
    write_csv(file = paste0(data_path,'map_iso_gene_vM34.csv'), col_names = FALSE)
}
if(TRUE){
  t0 <- Sys.time()
  if(TRUE){
    tmp2 <- parLapply(cl, sample_names, function(sample_name) {
      load(paste0(data_path,'mouse_pep_mapped.Rdata'))
      PepAbund <- 
        tibble(
          Y  = PepMapped[[sample_name]],
          EC = PepMapped$iso_ids
        ) %>%
        filter(Y > 0) %>%
        as.data.frame()
      suppressMessages(
        pep_data_loaded <- IsoBayes::load_data(
          path_to_peptides_psm = PepAbund,
          input_type = "other",
          abundance_type = 'intensities',
          PEP = FALSE
        )
      )
      suppressMessages(
        isoRes <- IsoBayes::inference(
          loaded_data = pep_data_loaded,
          map_iso_gene = paste0(data_path,'map_iso_gene_vM34.csv')
        )
      )
      
      write_rds(isoRes, file = paste0(output_path,sample_name,'_isoRes.rds'))
      
      return(NULL)
    })
    
    t1 <- Sys.time()
    time_cost <- difftime(t1, t0, units = "mins")
    print(time_cost)
  }
}

stopCluster(cl)

## send mails after finishing
receiver <- "chunxuhan007@gmail.com"
sender <- "18501059799@163.com"
emailSubject <- paste0('The quantification for J and A is done!')
emailBody <- "The R script has completed its execution. Please go and check the status and run the following samples."
send.mail(from = sender,
          to = receiver,
          subject = emailSubject,
          body = emailBody,
          smtp = list(host.name="smtp.163.com",
                      port=465, 
                      user.name=sender, 
                      passwd="ASGIPAQWXGMAONQY", 
                      ssl=TRUE),
          authenticate = TRUE,
          send = TRUE,
          encoding = "utf-8" 
)

## merge the results
file_paths <- list.files(output_path, pattern = "\\.rds$", full.names = TRUE)

process_file <- function(file_path) {
  data_list <- readRDS(file_path)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  file_name <- sub("_isoRes$", "", file_name)
  
  isoform_results <- data_list$isoform_results %>%
    select(Isoform, Abundance) %>%
    dplyr::rename(!!file_name := Abundance)
  
  gene_abundance <- data_list$gene_abundance %>%
    select(Gene, Abundance) %>%
    dplyr::rename(!!file_name := Abundance)
  
  list(isoform_results = isoform_results, gene_abundance = gene_abundance)
}


all_results <- lapply(file_paths, process_file)

isoform_mtx <- Reduce(function(x, y) {
  merge(x, y, by = "Isoform", all = TRUE)
}, lapply(all_results, `[[`, "isoform_results"))

# 合并gene_abundance
gene_mtx <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)
}, lapply(all_results, `[[`, "gene_abundance"))

rownames(isoform_mtx) <- isoform_mtx$Isoform
isoform_mtx$Isoform <- NULL

rownames(gene_mtx) <- gene_mtx$Gene
gene_mtx$Gene <- NULL

save(isoform_mtx,gene_mtx,file = paste0(work_path,'hema_iso_gene_mtx.Rdata'))
