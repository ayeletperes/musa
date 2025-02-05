if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  data.table,
  ggplot2,
  patchwork,
  ggrepel,
  extrafont,
  ggseqlogo,
  dplyr,
  extrafont,
  ComplexUpset,
  ggfun
)

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))


## load the contigs

load("data/rhesus_macaque_data/digger_files/new_IGH_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGH <- chain_df
load("data/rhesus_macaque_data/digger_files/new_IGL_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGL <- chain_df
load("data/rhesus_macaque_data/digger_files/new_IGK_digger_df_filter.rda")
setDT(chain_df)
digger_df_IGK <- chain_df
# Combine data frames
digger_df <- rbind(digger_df_IGH, digger_df_IGL, digger_df_IGK, fill=T)

new_subject_names <- fread("data/rhesus_macaque_data/new_subject_names.csv")
new_subject_names_dict <- setNames(new_subject_names$new_name, new_subject_names$sample)

digger_df[,new_subject:=new_subject_names_dict[subject]]

digger_df <- digger_df[!is.na(new_subject)]
# convert the current asc to the newest version.
MUSA <- fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_set_with_novel_2025-02-05_with_asc.csv")


allele_seed_novel_to_new_allele <- setNames(MUSA$new_allele, MUSA$allele_seed_novel)
digger_df[,new_allele:=allele_seed_novel_to_new_allele[asc]]

filtered_data <- MUSA[sample_count_genomic != 0]
filtered_data <- filtered_data[, .(new_allele, gene_type, samples_genomic, samples_AIRRseq)]
genomic_expanded <- filtered_data[, .(subject = unlist(strsplit(samples_genomic, ","))), by = .(new_allele, gene_type)]
genomic_expanded[, tag := "genomic"]
genomic_expanded[, gene := sub("\\*.*", "", new_allele)]
genomic_expanded <- genomic_expanded[!is.na(subject)]
airrseq_expanded <- filtered_data[, .(subject = unlist(strsplit(samples_AIRRseq, ","))), by = .(new_allele, gene_type)]
airrseq_expanded[, tag := "airrseq"]
airrseq_expanded[, gene := sub("\\*.*", "", new_allele)]
airrseq_expanded <- airrseq_expanded[!is.na(subject)]
combined_data <- rbind(genomic_expanded, airrseq_expanded, fill = TRUE)
combined_data[, subject := trimws(subject)]  # Trim white spaces
combined_data[, gene := sub("\\*.*", "", new_allele)]
summary_table <- combined_data[, .(unique_new_alleles_count = uniqueN(new_allele)), 
                               by = .(gene, gene_type, subject, tag)]
summary_table[tag=="genomic", unique_new_alleles_count:=1]
table_1 <- summary_table[order(tag == "airrseq", decreasing = TRUE)][
  !duplicated(paste(gene, subject)), ][
    , .(unique_subject_count = uniqueN(subject)), by = .(gene, gene_type, unique_new_alleles_count, tag)][
      ,source:=paste0(unique(tag),collapse = "_"),by=.(gene, gene_type)][
        grepl("airrseq",source)]#[(tag=="genomic" & unique_subject_count <106) | tag=="airrseq"]

# top each gene with a new category of unseen up to 106
gene_summary <- table_1[, .(total_subject_count = sum(unique_subject_count)), by = gene]
unseen_entries <- gene_summary[total_subject_count < 106, .(gene, unique_new_alleles_count = "unseen", 
                                                            unique_subject_count = 106 - total_subject_count)]
final_summary <- rbind(table_1, unseen_entries, fill = TRUE)
final_summary[is.na(tag), tag:="unseen"]
final_summary[,gene_type:=substr(gene,1,4)]
final_summary[,zygousity_state:=max(as.numeric(unique(unique_new_alleles_count[tag=="airrseq"]))), 
              by=.(gene, gene_type)]
final_summary[,zygousity_state_tag:=unique_new_alleles_count]
final_summary[unique_new_alleles_count=="unseen",zygousity_state_tag:=paste0("unseen-",zygousity_state)]
final_summary[tag=="genomic",zygousity_state_tag:="genomic"]
final_summary[,total_subject_count:=sum(unique_subject_count[tag=="airrseq"]), 
              by=.(gene, gene_type, zygousity_state)]

contigs_asc_all <- rbindlist(lapply(final_summary[grepl("V",gene_type) & zygousity_state>=4, unique(gene)], function(g){
  
  contigs_asc <- digger_df[grepl(g, new_allele),.(
    contigs_count = length(unique(contig)), 
    contigs = paste0(unique(contig),collapse = ","),
    alleles = paste0(unique(new_allele),collapse = ",")), by = new_subject]
  contigs_asc[,asc:=g]
  contigs_asc
}))

asc_check <- "IGHV3-GYPG"

# get the information for the selected asc
contigs_asc <- digger_df[grepl(asc_check, new_allele),.( 
  alleles = ifelse(unique(sense)=="-", 
                   paste0(rev(unique(new_allele)),collapse = ","),
                   paste0(unique(new_allele),collapse = ","))), by = .(new_subject, contig)]

contigs_asc2 <- contigs_asc[,.(allele = paste0(unique(alleles),";")), by = .(new_subject)]
contigs_asc2 <- contigs_asc2[,.(allele=unlist(strsplit(allele,";"))),by=.(new_subject)]

n_subjects <- length(unique(contigs_asc2$new_subject))
max_contigs <- contigs_asc2[,.(n_contigs=.N),by=new_subject][,max(n_contigs)]
mat <- matrix("", nrow = n_subjects, ncol = max_contigs,
              dimnames = list(unique(contigs_asc2$new_subject),paste0("contig_",1:max_contigs)))

for(i in 1:n_subjects){
  alleles <- contigs_asc2[new_subject==unique(contigs_asc2$new_subject)[i],allele]
  alleles_count <- stringi::stri_count(alleles,regex=",") + 1
  alleles_ <- setNames(alleles_count, alleles)
  
  mat[i,1:length(alleles)] <- names(sort(alleles_, decreasing = T))
}

mat_df <- as.data.table(mat)
mat_df$subject <- unique(contigs_asc$new_subject)
mat_df_melt <- melt.data.table(mat_df, id.vars = "subject", variable.name = "contig", value.name = "allele")

collapse_alleles <- function(allele) unique(paste0(sort(unlist(strsplit(unique(allele),","))), collapse = ","))

# rve17, move the values from contig 2 to contig 1, and sort
mat_df_melt[subject==new_subject_names_dict["RVe17"] & contig=="contig_2",contig:="contig_1"]
mat_df_melt[subject==new_subject_names_dict["RVe17"] & contig=="contig_3",contig:="contig_2"]
mat_df_melt[subject==new_subject_names_dict["RVe17"] & contig=="contig_1",allele:=collapse_alleles(allele), by = .(subject,contig)]

# R1r11, move the values from contig 2 to contig 1, and sort
mat_df_melt[subject==new_subject_names_dict["RQr11"] & contig=="contig_2",contig:="contig_1"]
mat_df_melt[subject==new_subject_names_dict["RQr11"] & contig=="contig_3",contig:="contig_2"]
mat_df_melt[subject==new_subject_names_dict["RQr11"] & contig=="contig_1",allele:=collapse_alleles(allele), by = .(subject,contig)]

# K410, move the values from contig 2 to contig 1, and sort
mat_df_melt[subject==new_subject_names_dict["K410"] & contig=="contig_2",contig:="contig_1"]
mat_df_melt[subject==new_subject_names_dict["K410"] & contig=="contig_3",contig:="contig_1"]
mat_df_melt[subject==new_subject_names_dict["K410"] & contig=="contig_1",allele:=collapse_alleles(allele), by = .(subject,contig)]

mat_df_melt <- mat_df_melt[!duplicated(paste(subject, contig, allele)),]

# order the rows
#subjects <- mat_df_melt[contig=="contig_1",][order(allele),(subject)]

## now split the allele into multiple rows
mat_df_melt <- mat_df_melt[,.(allele=sort(unique(unlist(strsplit(allele,","))))),by=.(subject,contig)]
mat_df_melt[,x:=1:.N,by=.(subject,contig)]
mat_df_melt[,contig_count:=length(unique(contig)),by=subject]

allele_counts <- mat_df_melt[, .(
  contig_1_alleles = paste0(sort(allele[contig == "contig_1"]), collapse = ","),
  contig_2_count = sum(contig == "contig_2"),
  contig_2_alleles = paste0(sort(allele[contig == "contig_2"]), collapse = ",")
), by = subject]
subjects <- allele_counts[order(contig_1_alleles, contig_2_count, contig_2_alleles), subject]
mat_df_melt$subject <- factor(mat_df_melt$subject, levels = subjects)

labs_ <- c("Contig A", "Contig B")
names(labs_) <- c("contig_1", "contig_2")

p_tiles <- ggplot(mat_df_melt, 
                  aes(x = x, y = subject, fill = allele)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  facet_grid(~contig, labeller = ggplot2::labeller(contig = labs_)) +
  theme_classic(base_size = 26) +
  labs(y = "Subject",x="",fill = "Allele") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_fill_manual(values = hcl.colors(uniqueN(mat_df_melt$allele), palette = "Temps")) +
  guides(fill = guide_legend(ncol = 2))
  
p_tiles_legend <- get_legend(p_tiles)

p_tiles_final <- wrap_elements((p_tiles + guides(fill="none"))/ p_tiles_legend + plot_layout(heights = c(4, 0.5)))


values = list(
  "1" = "#F8766D",
  "2" = "#C49A00","3" = "#53B400",
  "4" = "#00C094","5" = "#00B6EB",
  "6" = "#A58AFF","7" = "#FB61D7",
  "genomic" = "white",
  "unseen-1"= "#E6E6E6",
  "unseen-2"= "#D3D3D3",
  "unseen-3"= "#BFBFBF",
  "unseen-4"= "#A8A8A8",
  "unseen-5"= "#8C8C8C",
  "unseen-6"= "#4D4D4D",
  "unseen-7"= "#262626")


final_summary[zygousity_state_tag=="genomic", zygousity_state_tag:="Genomic Only"]

values = c(
  setNames((hcl.colors(7, palette = "Dark3")),as.character(1:7)),
  "Genomic Only" = "white",
  setNames((hcl.colors(15, palette = "Grays", rev = T))[5:11],paste0("unseen-",1:7)))


final_summary$zygousity_state_tag <- factor(final_summary$zygousity_state_tag, levels = rev(names(values)))

facet_labels <- function(value) {
  paste0("\u2264", value)  # Create parsed label format
}

zygousity_legend <- c()
plots <- list()


for(gene_type_ in c("IGHV","IGKV","IGLV")){
  final_summary_orderd_plot <- final_summary[gene_type==gene_type_][order(zygousity_state, -total_subject_count)][,gene:=factor(gene, levels = unique(gene))]
  ## create the box above each zygousity state tag
  
  p1 <- ggplot(final_summary_orderd_plot,aes(x = gene, 
                                             y = unique_subject_count, 
                                             fill = zygousity_state_tag)) +
    geom_bar(stat = "identity", width = 0.75) +
    scale_fill_manual(values = values, drop = FALSE, breaks = grep("unseen", names(values), value = T, invert = T)) +
    labs(y = "# Unique Subjects", fill = "Zygosity State") +
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) + 
    facet_grid(~zygousity_state, scales = "free_x", 
               space="free_x", labeller = as_labeller(facet_labels)) +
    theme_classic(base_size = 26) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.length.x=unit(0,'lines'),
          axis.line.x = element_line(colour = "gray10"),
          legend.key = element_rect(colour = "black"), 
          panel.spacing = unit(0, "lines"),
          strip.text = element_text(size = 13),
          legend.direction = "horizontal")
  
  ## grab the legend of IGKV
  if(gene_type_=="IGKV"){
    zygousity_legend <- get_legend(p1 + theme(panel.background = element_rect(colour = "transparent"),
                                              plot.background = element_rect(colour = "transparent"), 
                                              legend.background = element_rect(colour = "transparent"),
                                              legend.box.background = element_rect(colour = "transparent")))
  }
  p1 <- p1  +
      guides(fill = "none")
  
  table_2 <- combined_data[tag=="airrseq",.(unique_new_alleles_count = uniqueN(new_allele)), by=.(gene, gene_type)]
  zygous_state <- setNames(final_summary_orderd_plot$zygousity_state,final_summary_orderd_plot$gene)
  table_2_plot <- table_2[gene_type==gene_type_][,gene:=factor(gene, levels = levels(final_summary_orderd_plot$gene))]
  table_2_plot[,zygousity_state:=zygous_state[as.character(gene)]]
  
  table_2_plot[,gene_color := apply(table_2_plot,1, function(x){
    gene = x["gene"]
    col <- ifelse(gene=="IGHV3-GYPG", "darkslateblue", "black")
    glue::glue("<span style='color:{col}'>{gene}</span>")
  })]
  
  p2 <- ggplot(table_2_plot,aes(x = gene_color, y = unique_new_alleles_count)) +
    geom_bar(stat = "identity",fill = "black", width = 0.75) +
    theme_classic(base_size = 26) +
    facet_grid(~zygousity_state, scales = "free_x", space="free_x", labeller = as_labeller(facet_labels)) + 
    theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x=unit(0,'lines'),
          strip.text = element_blank(),
          strip.background = element_blank()
          ) +
    labs(y = "# Unique Alleles", x = gene_type_) +
    scale_y_reverse(expand = c(0, 0), breaks = seq(5, max(table_2_plot$unique_new_alleles_count), by = 5)) 
  
  plots[[gene_type_]] <- wrap_elements(p1/ plot_spacer() / p2 + plot_layout(heights = c(4, -0.72 ,4.5)))
}

layout <- "
AD
BD
CD
"
# legend

final_plot <- 
  (plots[[1]] + 
     inset_element(ggpubr::as_ggplot(zygousity_legend) + theme_void(), 0.3, -0.1, 0, 1, 
                   align_to = "panel", on_top = T, clip = T, ignore_tag = T)) + plots[[2]] + plots[[3]] + (free(p_tiles_final)) + 
  plot_layout(design = layout, widths = c(4,1)) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 46))


ggsave("figures/figure5.pdf", final_plot, width = 45, height = 30, dpi = 600, limitsize = FALSE, device = cairo_pdf)


