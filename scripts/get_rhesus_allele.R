library(dplyr)

##########################################################################################################################################################

# Collapse data for alleles since MUSA has separate rows for alleles with different RSS

musa <- read.csv("/data/MUSA_with_rss_and_leader_and_asc_2025-02-05.csv")
head (musa)
musa_in <- musa[musa$sample_count_genomic > 0 & musa$sample_count_AIRRseq > 0 & grepl("IGHV|IGKV|IGLV",musa$gene_type) & complete.cases(musa$sample_count_genomic, musa$sample_count_AIRRseq), ]

names(musa_in)

musa_in_filt <- musa_in[,c("new_allele","gene_type","seq","samples_genomic","samples_AIRRseq","kimdb","rhgldb.","imgt","vrc","guo")] %>%
  distinct() %>%
  group_by(new_allele, seq) %>%
  summarise(
    samples_genomic = paste(unique(samples_genomic), collapse = ", "),
    samples_AIRRseq = paste(unique(samples_AIRRseq), collapse = ", "),
    num_samples_genomic = length(unlist(strsplit(samples_genomic, ", "))),
    num_samples_AIRRseq = length(unlist(strsplit(samples_AIRRseq, ", "))),
    gene_type = paste(unique(gene_type), collapse = ", "),
    kimdb = paste(unique(kimdb), collapse = ", "),
    rhgldb. = paste(unique(rhgldb.), collapse = ", "),
    imgt = paste(unique(imgt), collapse = ", "),
    vrc = paste(unique(vrc), collapse = ", "),
    guo = paste(unique(guo), collapse = ", "),
    .groups = "drop"
  )

##########################################################################################################################################################

# Supplementary Table - Get the top 5 blast hits for each allele and add annotations from MUSA and bNAb class

germlines <- c("IGHV1-2","IGHV1-46","IGHV1-69","IGHV1-8","IGHV2-5","IGHV3-11","IGHV3-15","IGHV3-20","IGHV3-21","IGHV3-30","IGHV3-33","IGHV4-34","IGHV4-39","IGHV4-4","IGHV4-59",
               "IGKV1-13","IGKV1-33","IGKV1-39","IGKV1-5","IGKV1-6","IGKV1-9","IGKV2-24","IGKV2-28","IGKV2D-20","IGKV3-15","IGKV3-20","IGKV3-40",
               "IGLV1-47","IGLV1-51","IGLV2-11","IGLV2-14","IGLV2-23","IGLV2-8","IGLV3-1","IGLV3-19","IGLV3-21","IGLV3-25","IGLV7-43")

blast <- read.table("All.blastout", header = F) # All.blastout - concatenate output for IGHV, IGKV and IGLV
names(blast) <- c("human_allele","rhesus_allele","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs")
blast$human_gene = gsub("\\*.*","",blast$human_allele)
head(blast)

keep <- blast$human_gene %in% germlines

blast_germlines <- blast[keep,]

blast_top5 <- blast_germlines %>%
  group_by(human_allele) %>%
  slice_head(n = 5) %>%
  ungroup()

blast_top5_alias <- merge(blast_top5,musa_in_filt,by.x="rhesus_allele",by.y="new_allele",all.x=T)

bnab_class <- read.table("bnab_class.txt",header = T) 
bnab_class$Gene <- gsub("VH","IGHV",bnab_class$Gene)
bnab_class$Gene <- gsub("VK","IGKV",bnab_class$Gene)
bnab_class$Gene <- gsub("VL","IGLV",bnab_class$Gene)

bnab_class$Gene <- gsub("\u2013", "-", bnab_class$Gene)
blast_top5_alias$human_gene <- gsub("\u2013", "-", blast_top5_alias$human_gene)

bnab_class$Gene[!bnab_class$Gene %in% unique(blast_top5_alias$human_gene)]
# IGKV1-33 IGKV2D-20 IGKV3-40 no hits 

blast_top5_alias_class <- merge(blast_top5_alias, bnab_class, by.x = "human_gene", by.y = "Gene")

write.table(blast_top5_alias_class,"bnab_top5_annotated.txt", quote = F, sep="\t", row.names = F) 


###################################################################################################################################################

# Main table

blast_max <- blast_germlines %>%
  arrange(human_gene, desc(bitscore)) %>%  # Sort by human_gene and bitscore (descending)
  group_by(human_gene) %>%
  filter(bitscore == max(bitscore)) %>%  # Keep all rows with max bitscore
  ungroup()

blast_max_alias <- merge(blast_max,musa_in_filt,by.x="rhesus_allele",by.y="new_allele",all.x=T)
head(blast_max_alias)

blast_max_alias$human_gene <- gsub("\u2013", "-", blast_max_alias$human_gene)
blast_max_alias_class <- merge(blast_max_alias, bnab_class, by.x = "human_gene", by.y = "Gene")

write.table(blast_max_alias_class,"bnab_max_annotated.txt", quote = F, sep="\t", row.names = F) 
