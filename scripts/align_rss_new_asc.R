## align the short rss sequence.

if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  data.table,
  dplyr,
  tidyr,
  msa,
  Biostrings
)

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

# Compute the distance matrix between sequences.
# @param sequences A character vector of sequences.
# @param method The method to compute the distance. Either "lv" for Levenshtein distance or "hamming" for Hamming distance.
# @param trimm_prim_3 The number of nucleotides to trim from the 3' end of the sequences.
# @param type The type of output. Either "both" for both the dis and a matrix, "matrix" or "dist.
compute_reference_distance <- function(sequences, method = "lv", trimm_prim_3 = NULL, type = "both", quiet = FALSE){
  
  ## asset that trimm_prim_3 is either null or numeric
  if(!is.null(trimm_prim_3) && !is.numeric(trimm_prim_3)){
    stop("trimm_prim_3 must be either NULL or numeric")
  }
  ## assert that trimm_prim_3 is larger than 1
  if(!is.null(trimm_prim_3) && trimm_prim_3 < 1){
    stop("trimm_prim_3 must be larger than 1")
  }
  ## assert that the method is either lv or hamming
  if(!(method %in% c("lv", "hamming"))){
    stop("method must be either lv or hamming")
  }
  
  ## if trimm_prim_3 is not null, trimm the sequences
  if(!is.null(trimm_prim_3)){
    sequences <- substr(sequences, 1, trimm_prim_3)
  }
  
  ## if method is hamming make sure that all sequences have equal length. if not add N to max length and through a message
  if(method == "hamming"){
    if(!quiet) message(sprintf("sequences have been padded with N to the max length (%s).", max(nchar(sequences))))
    sequences <- gsub("\\s", "N", format(sequences, width = max(nchar(sequences))))
  }
  
  ## compute the distance matrix
  dist <- stringdist::stringdistmatrix(sequences, method = method, useNames = "names")
  
  if(type == "both"){
    return(list(dist = dist, matrix = as.matrix(dist)))
  }else{
    if(type == "matrix"){
      return(as.matrix(dist))
    }else{
      return(dist)
    }
  }
}

align_target_to_reference <- function(target, reference) {
  target_split <- strsplit(target, "")[[1]]
  reference_split <- strsplit(reference, "")[[1]]
  
  # Ensure both sequences are the same length
  n <- max(length(target_split), length(reference_split))
  target_split <- c(target_split, rep("-", n - length(target_split)))
  reference_split <- c(reference_split, rep("-", n - length(reference_split)))
  
  # Compute mismatches without gaps
  mismatches <- sum(target_split != reference_split)
  
  # Try introducing a single gap at all possible positions
  best_alignment <- target_split
  best_mismatches <- mismatches
  
  for (gap_pos in 1:n) {
    # Insert gap in target
    temp_target <- c(target_split[1:(gap_pos - 1)], "-", target_split[gap_pos:n])
    temp_target <- temp_target[1:n]  # Ensure length remains the same
    
    # Compute mismatches with this gap
    temp_mismatches <- sum(temp_target != reference_split)
    
    # Update best alignment if fewer mismatches
    if (temp_mismatches < best_mismatches) {
      best_alignment <- temp_target
      best_mismatches <- temp_mismatches
    }
  }
  
  # Return aligned target and reference
  list(
    target = paste(best_alignment, collapse = ""),
    reference = paste(reference_split, collapse = "")
  )
}

align_target_to_references <- function(target, references) {
  target_split <- strsplit(target, "")[[1]]
  
  best_alignment <- NULL
  best_mismatches <- Inf
  best_reference <- NULL
  
  # Iterate over each reference sequence
  for (reference in references) {
    reference_split <- strsplit(reference, "")[[1]]
    
    # Ensure both sequences are the same length
    n <- max(length(target_split), length(reference_split))
    padded_target <- c(target_split, rep("-", n - length(target_split)))
    padded_reference <- c(reference_split, rep("-", n - length(reference_split)))
    
    # Compute mismatches without gaps
    mismatches <- sum(padded_target != padded_reference)
    
    # Try introducing a single gap at all possible positions
    best_alignment_for_ref <- padded_target
    best_mismatches_for_ref <- mismatches
    
    for (gap_pos in 1:n) {
      # Insert gap in target
      temp_target <- c(padded_target[1:(gap_pos - 1)], "-", padded_target[gap_pos:n])
      temp_target <- temp_target[1:n]  # Ensure length remains the same
      
      # Compute mismatches with this gap
      temp_mismatches <- sum(temp_target != padded_reference)
      
      # Update best alignment if fewer mismatches
      if (temp_mismatches < best_mismatches_for_ref) {
        best_alignment_for_ref <- temp_target
        best_mismatches_for_ref <- temp_mismatches
      }
    }
    print(best_alignment_for_ref)
    # Update overall best alignment if this reference is better
    if (best_mismatches_for_ref < best_mismatches) {
      best_alignment <- best_alignment_for_ref
      best_mismatches <- best_mismatches_for_ref
      best_reference <- reference_split
    }
  }
  
  # Return the best alignment and reference
  list(
    target = paste(best_alignment, collapse = ""),
    reference = paste(best_reference, collapse = "")
  )
}

extract_rss <- function(data, cols, spacer_size = c(23,22), file_prefix = "tmp"){
  canonical_size <- intersect(c(23,12), spacer_size)
  rss_sequences_data <- data[!grepl("NONAMER not found|RSS not found", digger_notes, ignore.case = TRUE) & 
                               sample_count_genomic>0 & sample_count_AIRRseq>0 & sure_subject_count>0, ]
  rss_sequences_data[, rss := do.call(paste, c(.SD, sep = ",")), .SDcols = cols]
  
  rss_sequences_data[,gene_type := substr(new_allele, 1, 4)]
  rss_sequences_data <- rss_sequences_data[,.(new_allele, iglabel, gene_type, rss, sure_subject, sure_subject_count, l_part1, l_part2)]
  
  rss_sequences_data <- melt.data.table(rss_sequences_data, 
                                         id.vars = c("new_allele", "iglabel","gene_type","sure_subject","sure_subject_count","l_part1", "l_part2"))
  
  names(rss_sequences_data) <- c("new_allele", "iglabel","gene_type","samples","sure_sample_count", "l_part1", "l_part2", "rss_type", "rss_sequence")
  setDT(rss_sequences_data)
  
  rss_size <- spacer_size + 7 + 9 + 2
  rss_sequences_data <- rss_sequences_data[rss_sequence!="," & nchar(rss_sequence) %in% rss_size,.(new_allele, iglabel, gene_type, rss_type, rss_sequence, samples,sure_sample_count, l_part1, l_part2)]
  
  rss_sequences_data[, c("Heptamer", "Spacer", "Nonamer") := tstrsplit(rss_sequence, ",", fixed=TRUE)]
  rss_sequences_data <- unique(rss_sequences_data)
  rss_sequences_data[,rss:=gsub(",","",rss_sequence)]
  rss_sequences_data[,asc:=alakazam::getGene(new_allele, strip_d = F, omit_nl = F)]
  rss_sequences_data[,asc_tag:=gsub("IG[HKL]","",asc)]
  rss_sequences_data[,family:=unlist(lapply(strsplit(asc,"-"),'[[',1))]

  if(all(nchar(rss_sequences_data$Spacer)==canonical_size)){
    rss_sequences_data[,rss_aligned:=rss]
    return(rss_sequences_data)
  }
  
  spacer_sequences <- setNames(unique(rss_sequences_data$Spacer),unique(rss_sequences_data$Spacer))
  spacer_sequences_matrix <- compute_reference_distance(spacer_sequences, method = "lv", trimm_prim_3 = NULL, type = "both")
  id_of_short <- which(nchar(spacer_sequences)<canonical_size)
  
  ## group the non canonical rss by family and get the most common reference
  short_spacer_family <- rss_sequences_data[Spacer %in% names(id_of_short), .(Spacer, family)][!duplicated(Spacer)]
  
  references <- sapply(unique(short_spacer_family$family), function(fam){
    spacers <- short_spacer_family$Spacer[short_spacer_family$family==fam]
    range_ <- colnames(spacer_sequences_matrix$matrix)
    range_ <- range_[!range_ %in% names(id_of_short)]
    references <- unlist(lapply(spacers, function(i){
      x <- spacer_sequences_matrix$matrix[i,range_]
      x <- names(x)[x<=(min(x)+2)]
      x
    }))
    tab <- table(references)
    names(tab[tab==max(tab)])
  }, simplify = F)
  
  aligned_seq <- sapply(1:nrow(short_spacer_family), function(i){
    target_spacer <- short_spacer_family$Spacer[i]
    reference_spacers <- references[[short_spacer_family$family[i]]]
    alignment <- align_target_to_references(target = target_spacer, references = reference_spacers)
    alignment$target
  })
  names(aligned_seq) <- short_spacer_family$Spacer
  rss_sequences_data[,Spacer_aligned:=Spacer]
  rss_sequences_data[Spacer%in%names(aligned_seq),Spacer_aligned:=aligned_seq[Spacer]]
  rss_sequences_data[,rss_aligned:=paste0(Heptamer,Spacer_aligned,Nonamer)]

  return(rss_sequences_data)
}

macaque_alleles <- fread("data/MUSA_data/seed_set_with_novel_rss_filtered_2025-02-05_with_asc.csv")
allele_IGH <- macaque_alleles[chain == "IGH"]

## IGHV
allele_IGHV <- extract_rss(data = allele_IGH, cols = c("v_heptamer", "spacer_3", "v_nonamer"), 
                           spacer_size = c(23,22,21))
allele_IGHV[,rss_type:=paste0(rss_type,"_23")]

## IGHJ

allele_IGHJ <- extract_rss(data = allele_IGH, cols = c("j_heptamer", "spacer_5", "j_nonamer"), 
                           spacer_size = c(23,22))
allele_IGHJ[,rss_type:=paste0(rss_type,"_23")]

## IGHD
allele_IGHD_5 <- extract_rss(data = allele_IGH, cols = c("d_5_heptamer", "spacer_5", "d_5_nonamer"), 
                           spacer_size = c(12,11,13))
allele_IGHD_5[,gene_type:="IGHD_5"]
allele_IGHD_5[,rss_type:=paste0(rss_type,"_12")]

allele_IGHD_3 <- extract_rss(data = allele_IGH, cols = c("d_3_heptamer", "spacer_3", "d_3_nonamer"), 
                             spacer_size = c(12,11,13))
allele_IGHD_3[,gene_type:="IGHD_3"]
allele_IGHD_3[,rss_type:=paste0(rss_type,"_12")]

### IGK
allele_IGK <- macaque_alleles[chain == "IGK"]
allele_IGK_spacer <- melt.data.table(allele_IGK[sure_subject_count>0,.(gene_type, spacer_3, spacer_5)], id.vars = "gene_type")
allele_IGK_spacer[,len:=nchar(value)]

## IGKV

allele_IGKV <- extract_rss(allele_IGK, c("v_heptamer", "spacer_3", "v_nonamer"), 
                           spacer_size = c(12,11), "IGKV")
allele_IGKV[,rss_type:=paste0(rss_type,"_12")]

## IGKJ

allele_IGKJ <- extract_rss(allele_IGK, c("j_heptamer", "spacer_5", "j_nonamer"), 
                           spacer_size = c(23,22), "IGKJ")
allele_IGKJ[,rss_type:=paste0(rss_type,"_23")]

### IGL

allele_IGL <- macaque_alleles[chain == "IGL"]

allele_IGL_spacer <- melt.data.table(allele_IGL[sure_subject_count>0,.(gene_type, spacer_3, spacer_5)], id.vars = "gene_type")
allele_IGL_spacer[,len:=nchar(value)]

## IGKV

allele_IGLJ <- extract_rss(allele_IGL, c("j_heptamer", "spacer_5", "j_nonamer"), 
                           spacer_size = c(12,11), "IGLJ")
allele_IGLJ[,rss_type:=paste0(rss_type,"_12")]

## IGLV

allele_IGLV <- extract_rss(allele_IGL, c("v_heptamer", "spacer_3", "v_nonamer"), 
                           spacer_size = c(23,22), "IGLV")
allele_IGLV[,rss_type:=paste0(rss_type,"_23")]


## bind all tables

rss_sequences_data <- rbind(allele_IGHV, allele_IGHJ, allele_IGHD_5, allele_IGHD_3, 
                            allele_IGKV, allele_IGKJ, allele_IGLV, allele_IGLJ, fill=T)

fwrite(rss_sequences_data, file = sprintf("MUSA/seed_set_with_novel_rss_filtered_aligned_%s.tsv",format(Sys.time(), "%Y-%m-%d")), sep = "\t")

### sort human data.

process_rss_sequences <- function(data, cols, spacer_size, subject_count_filter, gene_tag, rss_tag, gapExtension = 500) {
  
  data_sub <- data[grepl(gene_tag, allele) & !grepl("NONAMER not found|RSS not found", notes, ignore.case = TRUE),]
  data_sub[, rss := do.call(paste, c(.SD, sep = ",")), .SDcols = cols]
  data_sub[, `:=`(
    sure_subject_count = ifelse(length(unique(rss)) == 1, 1, 0),
    non_sure_count = ifelse(length(unique(rss)) > 1, 1, 0),
    sure_subject = ifelse(length(unique(rss)) == 1, sample_name, NA)
  ), by = .(sample_name, allele)]
  
  rss_sequences_data <- data_sub[, .(
    gene_type = substr(allele, 1, 4),
    #samples = paste0(unique(sample_name), collapse = ","),
    #sample_count = length(unique(sample_name)),
    samples = paste0(unique(sure_subject[!is.na(sure_subject)]), collapse = ","),
    sample_count = sum(sure_subject_count[!duplicated(sample_name)]),
    non_sure_count = sum(non_sure_count[!duplicated(sample_name)])
  ), .(allele, rss)][sample_count > subject_count_filter, 
                     .(allele, gene_type, rss, samples, sample_count)]
  
  rss_sequences_data <- melt.data.table(
    rss_sequences_data, 
    id.vars = c("allele", "gene_type", "samples", "sample_count"),
    variable.name = "rss_type", value.name = "rss_sequence"
  )
  
  rss_size <- spacer_size + 7 + 9 + 2
  rss_sequences_data <- rss_sequences_data[
    rss_sequence != "," & nchar(rss_sequence) %in% rss_size, 
    .(allele, gene_type, rss_type, rss_sequence, samples, sample_count)
  ]
  
  rss_sequences_data[, c("Heptamer", "Spacer", "Nonamer") := tstrsplit(rss_sequence, ",", fixed = TRUE)]
  
  ## filter any heptamer or nonamer that are not 7 or 9
  rss_sequences_data <- rss_sequences_data[
    nchar(Heptamer) == 7 & nchar(Nonamer) == 9,
  ]
  
  rss_sequences_data <- unique(rss_sequences_data)
  rss_sequences_data[, rss := gsub(",", "", rss_sequence)]
  
  rss_sequences_data[, asc := alakazam::getGene(allele, strip_d = FALSE, omit_nl = FALSE)]
  if(grepl("J", gene_tag)){
    rss_sequences_data[, asc_tag := gsub(gene_tag, "", asc)]
  }else{
    rss_sequences_data[, asc_tag := unlist(lapply(strsplit(asc, "-"), '[[', 2))]  
  }
  
  
  rss_sequences_data[,rss_type:=paste0(rss_type,rss_tag)]
  
  rss_sequences_data[,rss_aligned:=rss]
  for(asc_ in unique(rss_sequences_data$asc)){
    if(all(nchar(rss_sequences_data[asc==asc_, Spacer])==max(spacer_size))){
      next
    }else{
      rss_sequences <- rss_sequences_data[asc==asc_, unique(rss)]
      alignment <- msa(DNAStringSet(rss_sequences), method = "ClustalW", verbose = F, gapExtension = gapExtension)
      rss_sequences_aligned <- setNames(as.character(alignment), rss_sequences)
      rss_sequences_data[asc==asc_,rss_aligned:=rss_sequences_aligned[rss]]
    }
  }
  
  return(rss_sequences_data)
}

df_human_igh <- fread("data/human_data/P25_imported_genes.csv")
setnames(df_human_igh, old = c(
  "vdjbase_allele","V-HEPTAMER", "V-SPACER", "V-NONAMER", "V-REGION-GAPPED", 
  "D-5_NONAMER", "D-5_SPACER", "D-5_HEPTAMER", "D-3_NONAMER", "D-3_SPACER", "D-3_HEPTAMER","D-REGION",
  "J-HEPTAMER", "J-SPACER", "J-NONAMER", "J-REGION"
), new = c(
  "allele","v_heptamer", "v_spacer","v_nonamer","v_sequence_gap",
  "d_5_nonamer", "d_5_spacer","d_5_heptamer","d_3_nonamer", "d_3_spacer","d_3_heptamer","d_sequence",
  "j_heptamer", "j_spacer","j_nonamer","j_sequence"
))

df_human_igh <- df_human_igh[,.(
  sample_name, allele, v_heptamer, v_spacer, v_nonamer, d_5_nonamer, d_5_spacer, d_5_heptamer,
  d_3_nonamer, d_3_spacer, d_3_heptamer, j_nonamer, j_spacer, j_heptamer, notes
)]

igh_ref <- c(
  names(tigger::readIgFasta("data/human_data/reference/IGHV.fasta")),
  names(tigger::readIgFasta("data/human_data/reference/IGHD.fasta")),
  names(tigger::readIgFasta("data/human_data/reference/IGHJ.fasta")))

df_human_ighv <- process_rss_sequences(df_human_igh[allele %in% igh_ref, ], c("v_heptamer", "v_spacer","v_nonamer"),
                                       gene_tag = "IGHV",spacer_size = c(23,22), subject_count_filter = 0, rss_tag = "_23")

df_human_ighd_5 <- process_rss_sequences(df_human_igh[allele %in% igh_ref, ], c("d_5_heptamer", "d_5_spacer","d_5_nonamer"),
                                       gene_tag = "IGHD",spacer_size = c(11,12), subject_count_filter = 0, rss_tag = "_12")
df_human_ighd_5[,gene_type:="IGHD_5"]

df_human_ighd_3 <- process_rss_sequences(df_human_igh[allele %in% igh_ref, ], c("d_3_heptamer", "d_3_spacer","d_3_nonamer"),
                                       gene_tag = "IGHD",spacer_size = c(11,12), subject_count_filter = 0, rss_tag = "_12")
df_human_ighd_3[,gene_type:="IGHD_3"]

df_human_ighj <- process_rss_sequences(data = df_human_igh[allele %in% igh_ref, ], cols = c("j_heptamer", "j_spacer","j_nonamer"),
                                       gene_tag = "IGHJ",spacer_size = c(23,22), subject_count_filter = 0, rss_tag = "_23")

df_human_igh_rss <- rbind(
  df_human_ighv,
  df_human_ighd_5,
  df_human_ighd_3,
  df_human_ighj
)

fwrite(df_human_igh_rss, file = "data/human_data/IGH_rss_sequences_data.tsv")

df_human_igk <- fread("data/human_data/PMID38844673_IGK_imported_genes.csv")
df_human_igk <- df_human_igk[df_human_igk$`Fully_Spanning_Reads_100%_Match`>2, ]
setnames(df_human_igk, old = c(
  "vdjbase_allele","V-HEPTAMER", "V-SPACER", "V-NONAMER", "V-REGION-GAPPED", 
  "D-5_NONAMER", "D-5_SPACER", "D-5_HEPTAMER", "D-3_NONAMER", "D-3_SPACER", "D-3_HEPTAMER","D-REGION",
  "J-HEPTAMER", "J-SPACER", "J-NONAMER", "J-REGION"
), new = c(
  "allele","v_heptamer", "v_spacer","v_nonamer","v_sequence_gap",
  "d_5_nonamer", "d_5_spacer","d_5_heptamer","d_3_nonamer", "d_3_spacer","d_3_heptamer","d_sequence",
  "j_heptamer", "j_spacer","j_nonamer","j_sequence"
))

igk_ref <- c(
  names(tigger::readIgFasta("data/human_data/reference/IGKV.fasta")),
  names(tigger::readIgFasta("data/human_data/reference/IGKJ.fasta")))


df_human_igkv <- process_rss_sequences(df_human_igk[allele %in% igk_ref, ], c("v_heptamer", "v_spacer","v_nonamer"),
                                       gene_tag = "IGKV",spacer_size = c(11,12), subject_count_filter = 0, rss_tag = "_12")

df_human_igkj <- process_rss_sequences(data = df_human_igk[allele %in% igk_ref, ], cols = c("j_heptamer", "j_spacer","j_nonamer"),
                                       gene_tag = "IGKJ",spacer_size = c(23,22), subject_count_filter = 0, rss_tag = "_23")


df_human_igl <- fread("data/human_data/P26_WGibson_213-sample-master_IGL.csv")
df_human_igl <- df_human_igl[df_human_igl$`Fully_Spanning_Reads_100%_Match`>2, ]
setnames(df_human_igl, old = c(
  "vdjbase_allele","V-HEPTAMER", "V-SPACER", "V-NONAMER", "V-REGION-GAPPED", 
  "D-5_NONAMER", "D-5_SPACER", "D-5_HEPTAMER", "D-3_NONAMER", "D-3_SPACER", "D-3_HEPTAMER","D-REGION",
  "J-HEPTAMER", "J-SPACER", "J-NONAMER", "J-REGION"
), new = c(
  "allele","v_heptamer", "v_spacer","v_nonamer","v_sequence_gap",
  "d_5_nonamer", "d_5_spacer","d_5_heptamer","d_3_nonamer", "d_3_spacer","d_3_heptamer","d_sequence",
  "j_heptamer", "j_spacer","j_nonamer","j_sequence"
))

igl_ref <- c(
  names(tigger::readIgFasta("data/human_data/reference/IGLV.fasta")),
  names(tigger::readIgFasta("data/human_data/reference/IGLJ.fasta")))


df_human_iglv <- process_rss_sequences(df_human_igl[allele %in% igl_ref, ], c("v_heptamer", "v_spacer","v_nonamer"),
                                       gene_tag = "IGLV",spacer_size = c(23,22), subject_count_filter = 0, rss_tag = "_23")

df_human_iglj <- process_rss_sequences(data = df_human_igl[allele %in% igl_ref, ], cols = c("j_heptamer", "j_spacer","j_nonamer"),
                                       gene_tag = "IGLJ",spacer_size = c(11,12), subject_count_filter = 0, rss_tag = "_12")



rss_sequences_data <- rbind(df_human_igh_rss, df_human_igkv, df_human_igkj, df_human_iglv, df_human_iglj)

fwrite(rss_sequences_data, file = "data/human_data/rss_sequences_data_aligned.tsv", sep = "\t")
