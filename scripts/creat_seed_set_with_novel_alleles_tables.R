# Load required libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(patchwork)
library(tidyr)
library(Biostrings)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(stringr)


# Load reference allele data from IGH, IGL, and IGK and combine into one dataframe
IGH <- read.csv("SeedSet/2025-01-03/IGH/IGH_crossref_with_digger_support.csv")
IGL <- read.csv("SeedSet/2025-01-03/IGL/IGL_crossref_with_digger_support.csv")
IGK <- read.csv("SeedSet/2025-01-03/IGK/IGK_crossref_with_digger_support.csv")

# Function to convert FASTA sequences into a dataframe
fasta_to_df <- function(fasta_file) {
  seq_names <- names(fasta_file)
  sequences <- as.character(fasta_file)
  return(data.frame(allele = seq_names, seq_gapped = sequences, stringsAsFactors = FALSE))
}

# Read the FASTA files (replace with actual paths if needed)
IGH_V_fasta <- readDNAStringSet("SeedSet/2025-01-03/IGH/V_gapped.fasta")
IGL_V_fasta <- readDNAStringSet("SeedSet/2025-01-03/IGL/V_gapped.fasta")
IGK_V_fasta <- readDNAStringSet("SeedSet/2025-01-03/IGK/V_gapped.fasta")

# Convert FASTA files to dataframes
IGH_V_df <- fasta_to_df(IGH_V_fasta)
IGL_V_df <- fasta_to_df(IGL_V_fasta)
IGK_V_df <- fasta_to_df(IGK_V_fasta)

# Merge the seq_gapped with the reference dataframes based on allele
IGH <- IGH %>% full_join(IGH_V_df, by = c("name" = "allele"))
IGL <- IGL %>% full_join(IGL_V_df, by = c("name" = "allele"))

IGK <- IGK %>% full_join(IGK_V_df, by = c("name" = "allele"))

reference_alleles <- bind_rows(IGH, IGL ,IGK)

reference_alleles <- reference_alleles %>%
  # Merge "rhgldb" and "cc_inferred" into a new column "rhgldb+"
  mutate(`rhgldb+` = case_when(
    rhgldb != "" ~ rhgldb,           # If rhgldb is not empty, use rhgldb
    cc_inferred != "" ~ cc_inferred, # If cc_inferred is not empty, use cc_inferred
    TRUE ~ ""             # Otherwise, set to NA
  )) %>%
  # Rename "kimdb.v1.1" to "kimdb"
  rename(`kimdb.v1.1` = "kimdb")%>%
  # Reorder the columns to place "rhgldb+" where "rhgldb" was
  select(name, gene_label, pubid, genbank, kimdb, `rhgldb+`, everything(), -rhgldb, -cc_inferred)

# Add a new column to mark entries present in the baseline reference
new_reference_df <- reference_alleles %>%
  dplyr::rename(
    allele = name)%>%
  mutate(
    in_baseline_reference = "TRUE"
  )

new_reference_df[sapply(new_reference_df, is.character)] <- 
  lapply(new_reference_df[sapply(new_reference_df, is.character)], 
         function(x) { ifelse(is.na(x) | x == "NA", "", x) })

selected_columns <- c(
  "asc", "asc_gene", "subject", "functional", "gene_type", "v_heptamer", "v_nonamer","j_heptamer", "j_nonamer", "d_3_heptamer" ,"d_3_nonamer","d_5_heptamer", "d_5_nonamer", "seq", "seq_gapped", "notes", "spacer_3", "spacer_5","l_part1","l_part2"
)

load("data/rhesus_macaque_data/digger_files/new_IGH_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGH <- chain_df[, ..selected_columns]
load("data/rhesus_macaque_data/digger_files/new_IGL_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGL <- chain_df[, ..selected_columns]
load("data/rhesus_macaque_data/digger_files/new_IGK_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGK <- chain_df[, ..selected_columns]
# Combine data frames

digger_df <- rbind(digger_df_IGH, digger_df_IGL, digger_df_IGK)


actual_df <- read.csv("data/rhesus_macaque_data/repertoire_genotype_all_04_02_rep_data.csv")




df_names <- read.csv("data/rhesus_macaque_data/new_subject_names.csv")

# Change `digger_df$subject == df_names$sample` to compare with `df_names$new_name`
digger_df$subject <- df_names$new_name[match(digger_df$subject, df_names$sample)]
# Change `actual_df$sample == df_names$sample` to compare with `df_names$new_name`
actual_df$sample <- df_names$new_name[match(actual_df$sample, df_names$sample)]


#actual_df <- actual_df %>%semi_join(digger_df, by = c("sample" = "subject", "allele" = "asc"))

#actual_df <- unique(actual_df)

m <- digger_df[order(digger_df$subject, digger_df$asc, digger_df$functional == "ORF"), ]

m <- m %>% group_by(subject, asc) %>% ungroup()
m <- m[!duplicated(m[c("subject", "asc" ,"v_heptamer", "v_nonamer","j_heptamer", "j_nonamer", "d_3_heptamer" ,"d_3_nonamer","d_5_heptamer", "d_5_nonamer", "seq", "seq_gapped", "functional", "spacer_3", "spacer_5","l_part1","l_part2")]), ]

# Create a flag for samples with both ORF and Functional alleles
m_with_flag <- m %>%
  group_by(subject, asc) %>%
  mutate(
    has_both = any(functional == "ORF") & any(functional == "Functional")
  ) %>%
  ungroup()

# Identify 'sure' and 'non-sure' subjects based on distinct values in heptamer and nonamer columns
m_with_subject_columns <- m_with_flag %>%
  group_by(subject, asc) %>%
  mutate(
    sure_subject = ifelse(
      !has_both & 
        n_distinct(v_nonamer) == 1 & n_distinct(v_heptamer) == 1 & 
        n_distinct(j_nonamer) == 1 & n_distinct(j_heptamer) == 1 & 
        n_distinct(d_3_heptamer) == 1 & n_distinct(d_3_nonamer) == 1 & 
        n_distinct(d_5_heptamer) == 1 & n_distinct(d_5_nonamer) == 1 & 
        n_distinct(spacer_3) == 1 & n_distinct(spacer_5) == 1 &
        n_distinct(l_part1) == 1 & n_distinct(l_part2) == 1, 
      subject, ""),
    non_sure_subject = ifelse(
      has_both | 
        n_distinct(v_nonamer) > 1 | n_distinct(v_heptamer) > 1 | 
        n_distinct(j_nonamer) > 1 | n_distinct(j_heptamer) > 1 | 
        n_distinct(d_3_heptamer) > 1 | n_distinct(d_3_nonamer) > 1 | 
        n_distinct(d_5_heptamer) > 1 | n_distinct(d_5_nonamer) > 1 | 
        n_distinct(spacer_3) > 1 | n_distinct(spacer_5) > 1 |
        n_distinct(l_part1) == 1 | n_distinct(l_part2) == 1, 
      subject, "")
  ) %>%
  ungroup()


actual_df$in_AIRRseq <- "TRUE"

m_with_subject_columns$in_genoimc <- "TRUE"


# Merge the actual dataframe with the subject columns from m_with_subject_columns
merged_df <- merge(actual_df, m_with_subject_columns, 
                   by.x = c("sample", "allele"), 
                   by.y = c("subject", "asc"), 
                   all.x = TRUE, all.y = TRUE)



# Select relevant columns for further analysis
allele_heptamer_nonamer <- merged_df[, c("sample", "allele","gene_type", 
                                         "v_heptamer", "v_nonamer", "j_heptamer", "j_nonamer",
                                         "d_3_heptamer" ,"d_3_nonamer","d_5_heptamer", "d_5_nonamer",
                                         "spacer_3", "spacer_5",
                                         "seq", "seq_gapped", "functional", 
                                         "sure_subject", "non_sure_subject", "notes", "z_score", "in_genoimc", "in_AIRRseq","l_part1","l_part2")]

allele_heptamer_nonamer[sapply(allele_heptamer_nonamer, is.character)] <- 
  lapply(allele_heptamer_nonamer[sapply(allele_heptamer_nonamer, is.character)], 
         function(x) { ifelse(is.na(x), "", x) })

# Summarize the data by allele and other characteristics
allele_heptamer_nonamer_count <- allele_heptamer_nonamer %>%
  group_by(across(-c(sample, notes, z_score, sure_subject,non_sure_subject, in_genoimc, in_AIRRseq))) %>%
  summarise(
    sample_count_genomic = n_distinct(sample[in_genoimc == "TRUE"]),
    samples_genomic = paste(unique(sample[in_genoimc == "TRUE"]), collapse = ", "),
    sample_count_AIRRseq = n_distinct(sample[in_AIRRseq == "TRUE"]),
    samples_AIRRseq = paste(unique(sample[in_AIRRseq == "TRUE"]), collapse = ", "),
    sure_subject_count = n_distinct(sure_subject[sure_subject != "" & in_AIRRseq == "TRUE"]),
    sure_subject = paste(unique(sure_subject[sure_subject != "" & in_AIRRseq == "TRUE"]), collapse = ", "),
    non_sure_subject_count = n_distinct(non_sure_subject[non_sure_subject != "" & in_AIRRseq == "TRUE"]),
    non_sure_subject = paste(unique(non_sure_subject[non_sure_subject != "" & in_AIRRseq == "TRUE"]), collapse = ", "),
    digger_notes = paste(unique(notes), collapse = "; "),
    z_score = paste(unique(z_score), collapse = "; "),
    .groups = 'drop'
  )

combined_with_references <- merge(allele_heptamer_nonamer_count,new_reference_df,by = c("allele", "seq", "seq_gapped"),all.x = TRUE, all.y = TRUE)

# Replace NA values in in_baseline_reference column with "FALSE"
combined_with_references$in_baseline_reference[is.na(combined_with_references$in_baseline_reference)] <- "FALSE"

# Replace NA in numeric columns with 0 and in character columns with empty strings
combined_with_references[sapply(combined_with_references, is.numeric)] <- 
  lapply(combined_with_references[sapply(combined_with_references, is.numeric)], 
         function(x) { ifelse(is.na(x) | x == "NA", 0, x) })

combined_with_references[sapply(combined_with_references, is.character)] <- 
  lapply(combined_with_references[sapply(combined_with_references, is.character)], 
         function(x) { ifelse(is.na(x) | x == "NA", "", x) })


# Add gene_type information by extracting the first four characters of the allele
combined_with_references <- combined_with_references %>% 
  mutate(asc = str_extract(allele, "(?<=-)[^*]+"),
         gene_type = substr(allele, 1, 4),
         chain = substr(allele, 1, 3)) %>%
  relocate(asc, gene_type, chain, .after = allele)

# Add RSS notes based on spacer_3 and spacer_5 sizes
combined_with_references <- combined_with_references %>%
  mutate(
    spacer_3_size = nchar(spacer_3),  # Calculate size of spacer_3
    spacer_5_size = nchar(spacer_5),  # Calculate size of spacer_5
    rss_notes = case_when(
      !(spacer_3_size %in% c(0, 12, 23)) & !(spacer_5_size %in% c(0, 12, 23)) ~ "Both spacer_3 and spacer_5 have variation at the size",
      !(spacer_3_size %in% c(0, 12, 23)) ~ "spacer_3 has variation at the size",
      !(spacer_5_size %in% c(0, 12, 23)) ~ "spacer_5 has variation at the size",
      TRUE ~ ""
    )
  )


#delet the incorect D.
median_value_df <- read.csv("data/rhesus_macaque_data/04_02_d_allele_more_0_median_value.csv")

# Filter alleles with median_value > 0.5
valid_alleles <- median_value_df %>%
  filter(median_value > 0.5) %>%
  pull(allele)  # Extract the list of valid alleles

# Remove incorrect rows
#combined_with_references <- combined_with_references %>%filter(!(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles))

combined_with_references <- combined_with_references %>%
  mutate(
    sample_count_AIRRseq = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, 0, sample_count_AIRRseq),
    sure_subject_count = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, 0, sure_subject_count),
    non_sure_subject_count = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, 0, non_sure_subject_count),
    samples_AIRRseq = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, "", samples_AIRRseq),
    sure_subject = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, "", sure_subject),
    non_sure_subject = ifelse(gene_type == "IGHD" & grepl("_", allele) & !allele %in% valid_alleles, "", non_sure_subject)
  )

# Update seq_gapped if allele contains "V", otherwise update seq
combined_with_references$seq_gapped <- ifelse(
  grepl("IGHV|IGLV|IGKV", combined_with_references$gene_type) & combined_with_references$seq_gapped == "",
  actual_df$ref[match(combined_with_references$allele, actual_df$allele)],
  combined_with_references$seq_gapped
)

combined_with_references$seq <-  ifelse(
  grepl("IGHV|IGLV|IGKV", combined_with_references$gene_type) & combined_with_references$seq == "",
  gsub("\\.", "", actual_df$ref[match(combined_with_references$allele, actual_df$allele)]),
  combined_with_references$seq
)

combined_with_references$seq <- ifelse(
  !grepl("V", combined_with_references$gene_type) & combined_with_references$seq == "",
  actual_df$ref[match(combined_with_references$allele, actual_df$allele)],
  combined_with_references$seq
)


filter_allele <- function(df) {
  df %>%
    # Filter only rows where allele contains 'V'
    filter(grepl("IGHV|IGLV|IGKV", allele)) %>%
    # Group by seq_gapped to identify duplicates
    group_by(seq) %>%
    # If there are duplicates, remove the novel allele with `_` if a non-novel allele exists
    filter(
      !(n_distinct(allele) > 1 & grepl("_", allele) & any(!grepl("_", allele)))
    ) %>%
    filter(
      !(n_distinct(allele) > 1 & (sample_count_AIRRseq > 0 & !(sample_count_genomic>0)) & any(sample_count_genomic>0))
    ) %>%
    # Ungroup after processing
    ungroup()
}

# Apply the function to your dataframe
v <- filter_allele(combined_with_references)

notv <- combined_with_references %>%
  filter(!grepl("IGHV|IGLV|IGKV", allele))

combined_with_references_update <- rbind(v,notv)


# Save the filtered dataframes for each chain as separate CSV files
unique_chains <- unique(combined_with_references_update$chain)

for (chain in unique_chains) {
  chain_df <- combined_with_references_update %>%
    filter(chain == !!chain)
  
  file_name <- paste0("data/rhesus_macaque_data/seed_set_with_novel_alleles/04_02_",chain, "_seed_and_novel_alleles.csv")
  write.csv(chain_df, file_name, row.names = FALSE)
}
