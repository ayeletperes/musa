if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  tigger,
  data.table,
  pbapply,
)

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

## creating the final bag of alleles including the baseline, novel allele, full IMGT set, VRC alleles, and Guo et al.

## baseline reference with new alleles

## Make sure we get the distribution iglabel and have all the baseline alleles
ref_old <- grep("V_iglabel[.]fasta",grep("_iglabel.fasta",grep("2025-01-03",list.files("../macaque_asc/work","fasta", full.names = T, recursive = T),value=T), value=T),invert=T, value=T)
names(ref_old) <- paste0("macaque_",gsub("_gapped","",gsub("_iglabel","",paste0(basename(dirname(dirname(ref_old))),tools::file_path_sans_ext(basename(ref_old))))))
ref_old <- ref_old[(!grepl("IG[KL]D", names(ref_old)))]
ref_old <- sapply(ref_old, readIgFasta)

ref_distribution_old <- grep("V[.]fasta",list.files("../macaque_asc/distributed/2025-01-03","fasta", full.names = T, recursive = T),invert = T,value = T)
names(ref_distribution_old) <- paste0("macaque_",gsub("_gapped","",paste0(basename(dirname(ref_distribution_old)),tools::file_path_sans_ext(basename(ref_distribution_old)))))
ref_distribution_old <- ref_distribution_old[(!grepl("IG[KL]D", names(ref_distribution_old)))]
ref_distribution_old <- sapply(ref_distribution_old, readIgFasta)


distribution_df_old <- rbindlist(lapply(ref_distribution_old, function(x){
  data.table(allele=names(x),sequence=gsub("[.]","",x))
}))

baseline_df_old <- rbindlist(lapply(ref_old, function(x){
  data.table(allele=names(x),sequence=gsub("[.]","",x))
}))

conversion_table <- merge.data.table(
  baseline_df_old, distribution_df_old, by = "sequence", all = T, suffixes = c("_baseline","_distributed")
)

conversion_list <- setNames(conversion_table$allele_baseline, conversion_table$allele_distributed)

## read the latest set of novel alleles and baseline alleles.
seed_set_and_novel <- rbindlist(
  list(
    "IGH" = fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/04_02_IGH_seed_and_novel_alleles.csv"),
    "IGK" = fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/04_02_IGK_seed_and_novel_alleles.csv"),
    "IGL" = fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/04_02_IGL_seed_and_novel_alleles.csv")
  ), fill = T)
seed_set_and_novel[seq_gapped=="",seq_gapped:=seq]
seed_set_and_novel[seq=="",seq:=gsub("[.]","",seq_gapped)]

seed_set_and_novel_cleaned <- seed_set_and_novel[,.(
  sample_count_AIRRseq = sum(sample_count_AIRRseq),
  sample_count_genomic = sum(sample_count_genomic),
  kimdb = unique(kimdb),
  `rhgldb+` = unique(`rhgldb+`),
  imgt = unique(imgt),
  trios = unique(trios),
  notes = unique(notes),
  iglabel = conversion_list[allele],
  novel = grepl("_",allele)
), by=.(allele,asc,gene_type,chain,seq,seq_gapped)]

seed_set_and_novel_cleaned[!grepl("_",allele),allele_seed:=allele]

fwrite(seed_set_and_novel_cleaned,
       sprintf("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_and_novel_alleles_with_iglabel_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")


# expected current set size
table(seed_set_and_novel_cleaned$gene_type)
# IGHD IGHJ IGHV IGKJ IGKV IGLJ IGLV 
# 74   19 1563    8  892   18  615 

MUSA_set <- copy(seed_set_and_novel_cleaned)

##############################
###### Add IMGT ##############
##############################

files <- list.files("../macaque_asc/other_ref_sets/imgt/",".fasta",recursive = T,full.names = T)
gapped_v <- grep("_gapped", files, value = T)
j_and_d <- grep("IG[HKL][DJ]", files, value = T)
files <- c(gapped_v, j_and_d)

# read the sets
imgt_sequences <- unlist(lapply(files, readIgFasta))
imgt_sequences_df <- data.table(allele = names(imgt_sequences), seq = gsub("[.]","",imgt_sequences), seq_gapped = imgt_sequences)

imgt_sequences_df[
  ,`:=`(
gene_type = substr(allele,1,4), 
chain = substr(allele,1,3)
)
]
# First check if any of the alleles that are in the set are not marked as imgt allele. If so mark them
imgt_sequences_not_marked <- imgt_sequences_df[seq %in% MUSA_set$seq[MUSA_set$imgt==""],]
imgt_sequences_not_marked <- setNames(imgt_sequences_not_marked$allele, imgt_sequences_not_marked$seq)

MUSA_set[seq %in% names(imgt_sequences_not_marked), imgt:=imgt_sequences_not_marked[seq], by=seq]

imgt_sequences_missing <- imgt_sequences_df[!seq %in% MUSA_set$seq,]
# Second, check if any are subset.

# special case is when the allele is a super or sub sequence of novel. then we will add it seperatly
missing_idx <- c()
for(i in 1:nrow(imgt_sequences_missing)){
  a <- imgt_sequences_missing$allele[i]
  s <- imgt_sequences_missing$seq[i]
  g_type <- imgt_sequences_missing$gene_type[i]
  idx <- MUSA_set[gene_type==g_type & grepl(s, seq),which = TRUE]
  if(length(idx)!=0){
    MUSA_set$imgt[idx] <- a
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("imgt: superseq of %s",a)  
      }else{
        if(any(!grepl(sprintf("imgt: superseq of %s",a),MUSA_set$notes[id]))){
          MUSA_set$notes[id] <- paste0(MUSA_set$notes[id],";",sprintf("imgt: superseq of %s",a))
        }
      }
    }
  }else{
    missing_idx <- c(missing_idx, i)
  }
}

## check if any of the missing is a super-sequence of any of the existing.
imgt_sequences_missing <- imgt_sequences_missing[missing_idx,]

missing_idx <- c()
for(i in 1:nrow(imgt_sequences_missing)){
  a <- imgt_sequences_missing$allele[i]
  s <- imgt_sequences_missing$seq[i]
  g_type <- imgt_sequences_missing$gene_type[i]
  seqs <- MUSA_set[gene_type==g_type,seq]
  idx <- sapply(seqs,function(x) grepl(x,s))
  idx <- which(idx)
  if(length(idx)!=0){
    idx <- MUSA_set[seq %in% names(idx),which = TRUE]
    #print(i)
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("imgt: subseq of %s",a)  
        MUSA_set$imgt[id] <- a
      }else{
        if(any(!grepl(sprintf("imgt: subseq of %s",a),MUSA_set$notes[id]))){
          notes <- c(MUSA_set$notes[id], sprintf("imgt: subseq of %s",a))
          MUSA_set$notes[id] <- paste0(unique(notes), collapse = ";")
          MUSA_set$imgt[id] <- paste0(unique(MUSA_set$imgt[id], a), collapse = ";")
        }
      }
    }
  }else{
    missing_idx <- c(missing_idx, i)
  }
}


imgt_sequences_missing <- imgt_sequences_missing[missing_idx,]
# Third, keep only sequences that are not already in the baseline set
## 148 alleles should be left. 
## add the sequences to the MUSA_set
# add missing columns
imgt_sequences_missing[,`:=`(
  asc = "", 
  sample_count_AIRRseq = 0, 
  sample_count_genomic = 0,
  kimdb = "", 
  imgt = allele, 
  `rhgldb+` = "", 
  trios = "", 
  novel = FALSE,
  notes = "The baseline didn't include IMGT alleles unless they were supported by another source",
  iglabel = NA
), by=allele]

# added alleles
# IGHD IGHV IGKV IGLV 
# 2  108   28   10 

MUSA_set <- rbind(MUSA_set, imgt_sequences_missing, fill=T)

##############################
###### Add VRC ###############
##############################

vrc_sequences <- fread("../macaque_asc/other_ref_sets/vrc/vrc_new_alleles.csv")
vrc_sequences <- vrc_sequences[,.(`Gene name`,Sequence)]
names(vrc_sequences) <- c("allele","seq")
vrc_sequences[,`:=`(
  gene_type = substr(allele,1,4), 
  chain = substr(allele,1,3)
)]

# IGHV IGKV IGLV 
# 121  163  111 

# First check if any of the alleles that are in the set are not marked as imgt allele. If so mark them
vrc_sequences_not_marked <- vrc_sequences[seq %in% MUSA_set$seq,]
vrc_sequences_not_marked <- setNames(vrc_sequences_not_marked$allele, vrc_sequences_not_marked$seq)

MUSA_set[,vrc:=""]
MUSA_set[seq %in% names(vrc_sequences_not_marked), vrc:=vrc_sequences_not_marked[seq], by=seq]

vrc_sequences_missing <- vrc_sequences[!seq %in% MUSA_set$seq,]
# Second, check if any are subset.
missing_idx <- c()
for(i in 1:nrow(vrc_sequences_missing)){
  a <- vrc_sequences_missing$allele[i]
  s <- vrc_sequences_missing$seq[i]
  g_type <- vrc_sequences_missing$gene_type[i]
  idx <- MUSA_set[gene_type==g_type & grepl(s, seq),which = TRUE]
  if(length(idx)!=0){
    MUSA_set$vrc[idx] <- a
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("vrc: superseq of %s",a)  
      }else{
        if(any(!grepl(sprintf("vrc: superseq of %s",a),MUSA_set$notes[id]))){
          MUSA_set$notes[id] <- paste0(MUSA_set$notes[id],";",sprintf("vrc: superseq of %s",a))
        }
      }
    }
    
  }else{
    missing_idx <- c(missing_idx, i)
  }
}


## check if any of the missing is a super-sequence of any of the existing.
vrc_sequences_missing <- vrc_sequences_missing[missing_idx,]
missing_idx <- c()
for(i in 1:nrow(vrc_sequences_missing)){
  a <- vrc_sequences_missing$allele[i]
  s <- vrc_sequences_missing$seq[i]
  g_type <- vrc_sequences_missing$gene_type[i]
  seqs <- MUSA_set[gene_type==g_type,seq]
  idx <- sapply(seqs,function(x) grepl(x,s))
  idx <- which(idx)
  if(length(idx)!=0){
    idx <- MUSA_set[seq %in% names(idx),which = TRUE]
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("vrc: subseq of %s",a)  
        MUSA_set$vrc[id] <- a
      }else{
        if(any(!grepl(sprintf("vrc: subseq of %s",a),MUSA_set$notes[id]))){
          notes <- c(MUSA_set$notes[id], sprintf("vrc: subseq of %s",a))
          MUSA_set$notes[id] <- paste0(unique(notes), collapse = ";")
          MUSA_set$vrc[id] <- paste0(unique(MUSA_set$vrc[id], a), collapse = ";")
        }
      }
    }
  }else{
    missing_idx <- c(missing_idx, i)
  }
}


vrc_sequences_missing <- vrc_sequences_missing[missing_idx,]

# Third, keep only sequences that are not already in the baseline set
## 79 alleles should be left. 
# add missing columns
vrc_sequences_missing[,`:=`(
  asc = "",
  sample_count_AIRRseq = 0, 
  sample_count_genomic = 0, 
  kimdb = "", 
  imgt = "", 
  `rhgldb+` = "", 
  trios = "", 
  vrc = allele,
  novel = FALSE,
  notes = "New source",
  iglabel = NA
), by=allele]

## as the sequences not aligned see if there are any close sequences, and use them for the alignment.
## We need aligned sequences for the asc inference

for(i in 1:nrow(vrc_sequences_missing)){
  a <- vrc_sequences_missing$allele[i]
  s <- vrc_sequences_missing$seq[i]
  g_type <- vrc_sequences_missing$gene_type[i]
  ref_seq <- MUSA_set[gene_type==g_type & nchar(seq)>=nchar(s),]
  ref_dist <- piglet:::allele_diff_indices_parallel2(ref_seq$seq, rep(s, nrow(ref_seq)), return_count = T)
  d <- ref_dist[which.min(ref_dist)]
  if(d < 10){
    gapped_ref <- ref_seq$seq_gapped[which.min(ref_dist)]
    s_gapped <- piglet:::insert_gaps2_vec(gapped_ref, s)
    vrc_sequences_missing[i,seq_gapped:=s_gapped]
  }else{
    #print(sprintf("allele %s, index %s, closest dist %s, closest allele %s", a, i, d, ref_seq$allele[which.min(ref_dist)]))
    gapped_ref <- ref_seq$seq_gapped[which.min(ref_dist)]
    s_gapped <- piglet:::insert_gaps2_vec(gapped_ref, s)
    vrc_sequences_missing[i,seq_gapped:=s_gapped]
  } 
}

# added alleles
# IGHV IGKV IGLV 
# 38   28   11 

MUSA_set <- rbind(MUSA_set, vrc_sequences_missing, fill=T)

##############################
###### Add Guo et al. ########
##############################

guo_sequences <- readIgFasta("../macaque_asc/other_ref_sets/guo/W_serise_RMDB.fasta")
guo_sequences <- data.table(allele = names(guo_sequences), seq = guo_sequences)

guo_sequences[,`:=`(
  gene_type = substr(allele,1,4), 
  chain = substr(allele,1,3)
)]

# IGHV 
# 197 

# First check if any of the alleles that are in the set are not marked as imgt allele. If so mark them
guo_sequences_not_marked <- guo_sequences[seq %in% MUSA_set$seq,]
guo_sequences_not_marked <- setNames(guo_sequences_not_marked$allele, guo_sequences_not_marked$seq)

MUSA_set[,guo:=""]
MUSA_set[seq %in% names(guo_sequences_not_marked), guo:=guo_sequences_not_marked[seq], by=seq]

# Second, keep only sequences that are not already in the baseline set
guo_sequences_missing <- guo_sequences[!seq %in% MUSA_set$seq,]
# Second, check if any are subset.
missing_idx <- c()
for(i in 1:nrow(guo_sequences_missing)){
  a <- guo_sequences_missing$allele[i]
  s <- guo_sequences_missing$seq[i]
  g_type <- guo_sequences_missing$gene_type[i]
  idx <- MUSA_set[gene_type==g_type & grepl(s, seq),which = TRUE]
  if(length(idx)!=0){
    MUSA_set$guo[idx] <- a
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("guo: superseq of %s",a)  
      }else{
        if(any(!grepl(sprintf("guo: superseq of %s",a),MUSA_set$notes[id]))){
          MUSA_set$notes[id] <- paste0(MUSA_set$notes[id],";",sprintf("guo: superseq of %s",a))
        }
      }
    }
    
  }else{
    missing_idx <- c(missing_idx, i)
  }
}

guo_sequences_missing <- guo_sequences_missing[missing_idx,]
missing_idx <- c()
for(i in 1:nrow(guo_sequences_missing)){
  a <- guo_sequences_missing$allele[i]
  s <- guo_sequences_missing$seq[i]
  g_type <- guo_sequences_missing$gene_type[i]
  seqs <- MUSA_set[gene_type==g_type,seq]
  idx <- sapply(seqs,function(x) grepl(x,s))
  idx <- which(idx)
  if(length(idx)!=0){
    idx <- MUSA_set[seq %in% names(idx),which = TRUE]
    for(id in idx){
      if(MUSA_set$novel[id]==TRUE){
        next
      }
      
      if(MUSA_set$notes[id]==""){
        MUSA_set$notes[id] <- sprintf("guo: subseq of %s",a)  
        MUSA_set$guo[id] <- a
      }else{
        if(any(!grepl(sprintf("guo: subseq of %s",a),MUSA_set$notes[id]))){
          notes <- c(MUSA_set$notes[id], sprintf("guo: subseq of %s",a))
          MUSA_set$notes[id] <- paste0(unique(notes), collapse = ";")
          MUSA_set$guo[id] <- paste0(unique(MUSA_set$guo[id], a), collapse = ";")
        }
      }
    }
  }else{
    missing_idx <- c(missing_idx, i)
  }
}

guo_sequences_missing <- guo_sequences_missing[missing_idx,]

## 136 alleles should be left. 

guo_sequences_missing[,`:=`(
  asc = "",
  sample_count_AIRRseq = 0, 
  sample_count_genomic = 0, 
  kimdb = "", 
  imgt = "", 
  `rhgldb+` = "", 
  trios = "", 
  vrc = "",
  guo = allele,
  novel = FALSE,
  notes = "New source",
  iglabel = NA
), by=allele]

## as the sequences not aligned see if there are any close sequences, and use them for the alignment.
## We need aligned sequences for the asc inferernce

for(i in 1:nrow(guo_sequences_missing)){
  a <- guo_sequences_missing$allele[i]
  s <- guo_sequences_missing$seq[i]
  g_type <- guo_sequences_missing$gene_type[i]
  ref_seq <- MUSA_set[gene_type==g_type & nchar(seq)>=nchar(s),]
  ref_dist <- piglet:::allele_diff_indices_parallel2(ref_seq$seq, rep(s, nrow(ref_seq)), return_count = T)
  d <- ref_dist[which.min(ref_dist)]
  if(d < 10){
    gapped_ref <- ref_seq$seq_gapped[which.min(ref_dist)]
    s_gapped <- piglet:::insert_gaps2_vec(gapped_ref, s)
    guo_sequences_missing[i,seq_gapped:=s_gapped]
  }else{
    #print(sprintf("allele %s, index %s, closest dist %s, closest allele %s", a, i, d, ref_seq$allele[which.min(ref_dist)]))
    gapped_ref <- ref_seq$seq_gapped[which.min(ref_dist)]
    s_gapped <- piglet:::insert_gaps2_vec(gapped_ref, s)
    guo_sequences_missing[i,seq_gapped:=s_gapped]
  } 
}

MUSA_set <- rbind(MUSA_set, guo_sequences_missing, fill=T)

# IGHD IGHJ IGHV IGKJ IGKV IGLJ IGLV 
# 76   19 1857    8  948   18  637 
###############################

# in cases where a novel allele was seen in any of the sources.
# we will change the allele annotation to the source allele.
# priority is IMGT, VRC, Guo.
# we will change the novel call to FALSE, but retain the "allele_old_distribution" and add a column novel_notes
# we need to remove the asc annotation,

MUSA_set_backup <- copy(MUSA_set)

MUSA_set[,allele_seed_novel:=allele]

MUSA_set[novel==TRUE & (imgt!="" |  vrc!="" | guo!=""),`:=`(
  novel_notes = "Seen in external sources",
  novel = FALSE,
  allele = apply(.SD, 1, function(x){
    if(x["imgt"]!=""){
      x["imgt"]
    }else if(x["vrc"]!=""){
      x["vrc"]
    }else if(x["guo"]!=""){
      x["guo"]
    }else{
      x["allele"]
    }
  }),
  asc = ""
)]


seed_to_new_novel <- setNames(MUSA_set$allele, MUSA_set$allele_seed_novel)

seed_set_and_novel[,allele_seed_novel:=allele]
seed_set_and_novel[,allele_update:=seed_to_new_novel[allele]]

fwrite(seed_set_and_novel,sprintf("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_set_and_novel_%s_update.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")


# get iglabel for the missing sequences, save it by source
save(MUSA_set, file = sprintf("data/MUSA_data/MUSA_set_all_%s.rds",format(Sys.time(), "%Y-%m-%d")))
fwrite(MUSA_set,sprintf("data/MUSA_data/MUSA_set_all_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")

# bosinger_watson
fwrite(MUSA_set[is.na(iglabel) & novel==TRUE,],sprintf("data/MUSA_data/without_iglabel/MUSA_set_alleles_without_iglabel_bosinger_watson_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")

# imgt
fwrite(MUSA_set[is.na(iglabel) & imgt!="" & novel!=TRUE,],sprintf("data/MUSA_data/without_iglabel/MUSA_set_alleles_without_iglabel_imgt_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")

# vrc
fwrite(MUSA_set[is.na(iglabel) & vrc!="" & novel!=TRUE,],sprintf("data/MUSA_data/without_iglabel/MUSA_set_alleles_without_iglabel_vrc_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")

# guo
fwrite(MUSA_set[is.na(iglabel) & guo!="" & novel!=TRUE,],sprintf("data/MUSA_data/without_iglabel/MUSA_set_alleles_without_iglabel_guo_%s.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")


## Add the missing iglabels. 

novel_allele_iglabel <- rbindlist(lapply(grep(format(Sys.time(), "%Y-%m-%d"),list.files("data/MUSA_data/with_iglabel","updated_with_iglabel",full.names = T),value=T),fread),fill=T)
any(is.na(novel_allele_iglabel$iglabel))
any(novel_allele_iglabel$iglabel=="")

## check if the novel is perhapse in the database
iglabel_database <- rbindlist(lapply(list.files("data/MUSA_data/db","_with_trio",full.names = T),fread),fill=T)

allele_set <- setNames(novel_allele_iglabel$iglabel, novel_allele_iglabel$seq_gapped)
MUSA_set[,label:=gsub("[*]00","",gsub("IG[HKL][VDJ]0-","",iglabel))]
MUSA_set[seq_gapped %in% names(allele_set),label:=allele_set[seq_gapped]]
MUSA_set[,iglabel:=paste0(gene_type,"0-",label,"*00")]
MUSA_set[,present:=sample_count_AIRRseq>0]
## check that all alleles have an iglabel
any(is.na(MUSA_set$iglabel))
any(grepl("IG[HKL][VDJ]0-NA[*]00",MUSA_set$iglabel))
any(grepl("IG[HKL][VDJ]0-[*]00",MUSA_set$iglabel))
MUSA_set[is.na(allele_seed), allele_seed:=""]
# remove duplicated
MUSA_set_collapsed <- MUSA_set[,.(
  allele = unique(allele)[1],
  asc = unique(asc),
  sample_count_AIRRseq = sum(sample_count_AIRRseq),
  sample_count_genomic = sum(sample_count_genomic),
  in_baseline_reference = !unique(novel),
  kimdb = paste0(unique(kimdb),collapse = ";"),
  `rhgldb+` = paste0(unique(`rhgldb+`),collapse = ";"),
  imgt = paste0(unique(imgt),collapse = ";"),
  trios = paste0(unique(trios),collapse = ";"),
  vrc = paste0(unique(vrc),collapse = ";"),
  guo = paste0(unique(guo),collapse = ";"),
  notes = paste0(unique(notes),collapse = ";"),
  iglabel = paste0(unique(iglabel),collapse = ";"),
  allele_seed = paste0(unique(allele_seed),collapse = ";"),
  novel = unique(novel),
  novel_notes = paste0(unique(novel_notes),collapse = ";"),
  present = unique(present)
),by=.(seq,seq_gapped,gene_type,chain)]

## check for duplicated iglabel

MUSA_set_collapsed[duplicated(iglabel) | duplicated(iglabel, fromLast=T)]

## there are 4 duplicates IGLV21-KXCY*01 & IGLV21-6YH6*03_g112a_a220g_a224g, and IGLJ9-PIRL*01 & IGLJ9-PIRL*01_i0ttggg

## we will keep both alleles  IGLV21-KXCY*01 (TCTGTGCTGACTCAGCCACCCTCCCTCTCTGCATCCCTGGAAGCATTAGCCAGACTCACCTGCACCCTGAGCAGTGGCATCAGTGTTGGTGGAAAAATTGTATACTGGTACCAGCAGAAGCCAGGGAGCAATCCCCGGTATCTCCTGAGCTACTACTCAGAGTCAAGTAAGCACCAGGGCTCTGGAGTCCCCGGCCGCTTCTCTGGATCCAAAGATGCCTCAACTAACTCAGGGATTCTGCACGTCTCTGGGCTGCAGCCTGAGGATGAGGCTGACTATTATTGTAAGATATGGCATGACAGCATTAATGCTTA) 
## & IGLV21-6YH6*03_g112a_a220g_a224g but give the longer one a different iglabel as both are seen in the cohort.

## we will remove IGLJ9-PIRL*01 (TGTTCGGCGAGGGGACCAAGCTGACCATCCTAG) as it is a shorter version of IGLJ9-PIRL*01_i0ttggg (TTGGGTGTTCGGCGAGGGGACCAAGCTGACCATCCTAG) and is not seen in any subject.

MUSA_set_collapsed_filtered <- MUSA_set_collapsed[seq!="TGTTCGGCGAGGGGACCAAGCTGACCATCCTAG",]

## write file with iglabels
fwrite(MUSA_set_collapsed_filtered,sprintf("data/MUSA_data/MUSA_set_all_%s_with_iglabel.csv",format(Sys.time(), "%Y-%m-%d")), sep = ",")



