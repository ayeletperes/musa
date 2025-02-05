## conversion table between ditrbiution reference 2024-06-24 and 2025-01-03

if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  data.table,
  pbapply,
)

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

## Make sure we get the distribution iglabel and have all the baseline alleles
iglabel_2025_01_03 <- grep("V_iglabel[.]fasta",grep("_iglabel.fasta",grep("2025-01-03",list.files("../macaque_asc/work","fasta", full.names = T, recursive = T),value=T), value=T),invert=T, value=T)
names(iglabel_2025_01_03) <- paste0("macaque_",gsub("_gapped","",gsub("_iglabel","",paste0(basename(dirname(dirname(iglabel_2025_01_03))),tools::file_path_sans_ext(basename(iglabel_2025_01_03))))))
iglabel_2025_01_03 <- iglabel_2025_01_03[(!grepl("IG[KL]D", names(iglabel_2025_01_03)))]
iglabel_2025_01_03 <- sapply(iglabel_2025_01_03, readIgFasta)


asc_2025_01_03 <- grep("V[.]fasta",list.files("../macaque_asc/distributed/2025-01-03","fasta", full.names = T, recursive = T),invert = T,value = T)
names(asc_2025_01_03) <- paste0("macaque_",gsub("_gapped","",paste0(basename(dirname(asc_2025_01_03)),tools::file_path_sans_ext(basename(asc_2025_01_03)))))
asc_2025_01_03 <- asc_2025_01_03[(!grepl("IG[KL]D", names(asc_2025_01_03)))]
asc_2025_01_03 <- sapply(asc_2025_01_03, readIgFasta)


## we need to keep the latest reference the AIRR-seq data was run with.
asc_2024_06_24 <- grep("V[.]fasta",list.files("../macaque_asc/distributed/2024-06-24","fasta", full.names = T, recursive = T),invert = T,value = T)
names(asc_2024_06_24) <- paste0("macaque_",gsub("_gapped","",paste0(basename(dirname(asc_2024_06_24)),tools::file_path_sans_ext(basename(asc_2024_06_24)))))
asc_2024_06_24 <- asc_2024_06_24[(!grepl("IG[KL]D", names(asc_2024_06_24)))]
asc_2024_06_24 <- sapply(asc_2024_06_24, readIgFasta)

asc_2025_01_03_df <- rbindlist(lapply(asc_2025_01_03, function(x){
  data.table(allele=names(x),sequence=gsub("[.]","",x))
}))

asc_2024_06_24_df <- rbindlist(lapply(asc_2024_06_24, function(x){
  data.table(allele=names(x),sequence=gsub("[.]","",x))
}))

iglabel_2025_01_03_df <- rbindlist(lapply(iglabel_2025_01_03, function(x){
  data.table(allele=names(x),sequence=gsub("[.]","",x))
}))

conversion_table <- merge.data.table(
  iglabel_2025_01_03_df, asc_2025_01_03_df, by = "sequence", all = T, suffixes = c("_iglabel","_2025-01-03")
)

conversion_table_with_run <- merge.data.table(
  conversion_table, asc_2024_06_24_df, by = "sequence", all = T
)

names(conversion_table_with_run)[4] <- "allele_2024-06-24"

conversion_list <- setNames(conversion_table_with_run[["allele_iglabel"]], conversion_table_with_run[["allele_2024-06-24"]])
conversion_list_new_distribution <- setNames(conversion_table_with_run[["allele_2025-01-03"]], conversion_table_with_run[["allele_2024-06-24"]])
conversion_oppiste_list <- setNames(conversion_table_with_run[["allele_2025-01-03"]], conversion_table_with_run[["allele_2024-06-24"]])

fwrite(conversion_table_with_run, "SeedSet/conversion_table_2024-06-24_to_2025-01-03.csv", sep = ",")
save(conversion_list, conversion_list_new_distribution, file = "SeedSet/conversion_list_2024-06-24_to_2025-01-03.rds")
