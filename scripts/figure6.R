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
  ComplexUpset
)
font_import(prompt = F)
loadfonts()

## correct to include the spacer for the count and replace it after with space.
## add distnace matrix of the complemet reverse. HM distance. take the minimum. take all the 23 together and all the 12 together. with the spacer. and add NJ tree. 

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

macaque_data = fread("data/MUSA_data/seed_set_with_novel_rss_filtered_aligned_2025-02-04.tsv")
macaque_data <- macaque_data[sure_sample_count>0,]
human_data = fread("data/human_data/rss_sequences_data_aligned.tsv")
rss_sequences_data <- rbindlist(list(
  "macaque" = macaque_data, 
  "human" = human_data), use.names = T, idcol = "specie", fill = T)
## remove the spacer and add N instead

## first we do the count.
rss_sequences_data[, rss_sequence := paste0(Heptamer, "NNN", Nonamer)]

rss_sequences <- rss_sequences_data[,.(count=.N), by = .(rss_type, specie, rss_aligned, gene_type, Nonamer, Heptamer, Spacer, rss_sequence)]

rss_sequences[,color:=ifelse(.N>1, "both", "single"), by = .(Nonamer, Heptamer, gene_type)]

rss_sequences_macaque <- rss_sequences[specie == "macaque",]


## check if how many of the rss are in imgt

imgt_rss <- tigger::readIgFasta("data/rhesus_macaque_data/imgt_rss.txt")

v_rss <- macaque_data[gene_type=="IGHV", unique(rss_aligned)]

print(sprintf("%s RSS overlap",length(intersect(v_rss, imgt_rss))))

print(sprintf("%s RSS are new",length(v_rss[!v_rss %in% imgt_rss])))

print(sprintf("%s RSS are unique to imgt",length(imgt_rss[!imgt_rss %in% v_rss])))


## check for how many alleles the rss pairing is new

print(sprintf("%s RSS are new pairs",macaque_data[rss_aligned %in% v_rss[!v_rss %in% imgt_rss],length(unique(new_allele))]))


## do the same for the leader as well. 
imgt_rss_leader <- tigger::readIgFasta("data/rhesus_macaque_data/imgt_rss_leader.txt")
macaque_alleles <- fread("data/MUSA_data/seed_set_with_novel_rss_filtered_2025-02-05_with_asc.csv")
MUSA_alleles <- fread("data/MUSA_data/MUSA_set_all_2025-02-05_with_asc.csv")
allele_seq <- setNames(macaque_alleles$seq, macaque_alleles$new_allele)
allele_novel <- setNames(MUSA_alleles$novel, MUSA_alleles$new_allele)

v_rss_leader <- copy(macaque_data[gene_type=="IGHV", ])
v_rss_leader[,seq:=allele_seq[new_allele]]
v_rss_leader[,novel:=allele_novel[new_allele]]

v_rss_leader[,in_imgt:=apply(v_rss_leader,1, function(x){
  
  l_part1 <- x[["l_part1"]]
  
  l_part2 <- x[["l_part2"]]
  
  seq <- x[["seq"]]
  
  rss_aligned <- x[["rss_aligned"]]
  
  ## first the rss - sequence combo
  rss_seq <- paste0(seq, rss_aligned)
  imgt_seq <- grep(rss_seq, imgt_rss_leader, value = T)
  
  if(length(imgt_seq)!=0){
    # check l_part1
    imgt_seq_l_part1 <- grep(l_part1, imgt_seq, value = T)
    if(length(imgt_seq_l_part1)!=0){
      # check l_part2
      imgt_seq_l_part1_l_part2 <- grep(l_part2, imgt_seq_l_part1, value = T)
      if(length(imgt_seq_l_part1_l_part2)!=0){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
  
})]

table(v_rss_leader$in_imgt)

print(sprintf("%s RSS-Leader are new",sum(!v_rss_leader$in_imgt)))

print(sprintf("%s RSS-Leader are new",v_rss_leader[novel==FALSE,sum(!in_imgt)]))

####


rss_sequences_macaque[,underline:=ifelse(length(unique(gene_type))>1, "underline", "no underline"), by = .(Nonamer, Heptamer)]

rss_sequences_macaque[,underline_full:=ifelse(length(unique(gene_type))>1, "underline", "no underline"), by = .(rss_aligned)]

rss_sequences_macaque[,spacer_share:=ifelse(length(unique(gene_type))>1, "underline", "no underline"), by = .(Spacer)]


rss_sequences_macaque <- rss_sequences_macaque[order(gene_type, -count, rss_aligned)]

rss_sequences_macaque[,ig_tag:=substr(gene_type,4,4)]

rss_sequences_macaque <- rss_sequences_macaque[order(ig_tag, -count)]

rss_sequences_macaque$tag <- ""

rss_names <- c()
i_v <- 0
i_j <- 0
i_d <- 0
for(gene_type_ in c("V","J","D")){
  seqs <- rss_sequences_macaque[ig_tag == gene_type_, unique(rss_sequence)]
  for(seq in seqs){
    
    if(seq %in% rss_names){
      next
    }
    
    if(gene_type_ == "V"){
      i_v <- i_v + 1
      rss_names[seq] <- paste0(gene_type_,i_v)
    }else{
      if(gene_type_ == "J"){
        i_j <- i_j + 1
        rss_names[seq] <- paste0(gene_type_,i_j)
      }else{
        i_d <- i_d + 1
        rss_names[seq] <- paste0(gene_type_,i_d)
      }
    }
    
  }
  
}


rss_sequences_macaque$tag <- rss_names[rss_sequences_macaque$rss_sequence]
rss_sequences_macaque[,rss_tag:=paste0(tag,"-",1:.N),by=tag]


rss_sequences_macaque[,nseq:=.N>15, by = .(gene_type)]
rss_sequences_macaque_subset <- rss_sequences_macaque[, head(.SD, 15),by=gene_type]
rss_sequences_macaque_subset <- rss_sequences_macaque_subset[order(gene_type, count, rss_sequence)]
rss_sequences_macaque_subset[, order := seq(1, by = 1, length.out = .N), by = .(gene_type)]
rss_sequences_macaque_subset[, count_char:= as.character(count)]

for(gene_type_ in c("IGHV","IGLV","IGKV","IGHD_3","IGHJ","IGLJ","IGKJ","IGHD_5")){
  
  tmp <- rss_sequences_macaque_subset[gene_type == gene_type_]
  
  if (tmp$nseq[1]) {
    new_rows <- tmp[rep(1, 3)]
    
    # Modify the new rows
    new_rows[, `:=`(
      order = 1:.N,
      rss_sequence = gsub("[ATCG]", " ", rss_sequence) %>%
        sub("N", " ", .) %>%
        sub("N", ".", .) %>%
        sub("N", " ", .),
      count_char = ".",
      underline = "no underline",
      color = "single",
      rss_tag = paste0(".",1:.N)
    )]
    
    rss_sequences_macaque_subset[gene_type == gene_type_, order := order + 3]
    rss_sequences_macaque_subset <- rbind(rss_sequences_macaque_subset, new_rows)
  }
  
}


rss_sequences_letters <- rss_sequences_macaque_subset[, .(letter = unlist(strsplit(rss_sequence, ""))), 
                                       by = .(gene_type, specie, rss_sequence, order, color, underline, rss_tag)]
rss_sequences_letters[,position := 1:(.N) , by = .(gene_type, rss_tag)]

rss_sequences_count <- rss_sequences_macaque_subset[, .(letter = count_char), 
                                     by = .(gene_type, specie, rss_sequence, order, color, underline, rss_tag)]
rss_sequences_count[,position := 1]

rss_sequences_tag <- rss_sequences_macaque_subset[, .(letter = rss_tag), 
                                             by = .(gene_type, specie, rss_sequence, order, color, underline, rss_tag)]
rss_sequences_tag[,position := 1]


scale_color_vector <- c(
  "single" = "#000000",
  "both" = "#000000",
  "underline" = "#000000"
)


scale_color_vector <- c(
  "single" = "#FF0000",
  "both" = "#000000",
  "underline" = "#000000"
)




max_height_larger_group <- max(rss_sequences_letters$order)
larger_group_id <- rss_sequences_letters$gene_type[which.max(rss_sequences_letters$order)]

diff_value_smaller_group <- sapply(c("IGHV","IGLV","IGKV","IGHD_3","IGHJ","IGLJ","IGKJ","IGHD_5"), function(seg){
  max_height_larger_group - max(rss_sequences_letters[gene_type == seg, order])
})

rss_by_gene_plot_list <- list()

for(gene_type_ in c("IGHV","IGLV","IGKV","IGHD_3","IGHJ","IGLJ","IGKJ","IGHD_5")){
  
  rss_sequences_data_segment <- rss_sequences_data[gene_type == gene_type_ & specie == "macaque",]
  
  if(gene_type_ %in% c("IGHD_3","IGHD_5","IGKV","IGLJ")){
    lab <- "12"
    if(any(nchar(rss_sequences_data_segment$Spacer)<12)){
      lab <- "11/\n12"
    }
  }else{
    lab <- "23"
    if(any(nchar(rss_sequences_data_segment$Spacer)<23)){
      lab <- "22/\n23"
    }
  }
  
  rss_sequences_letters_segment = rss_sequences_letters[gene_type == gene_type_,]
  rss_sequences_count_segment = rss_sequences_count[gene_type == gene_type_,]
  rss_sequences_tag_segment = rss_sequences_tag[gene_type == gene_type_,]
  
  diff <-diff_value_smaller_group[gene_type_]
  
  rss_sequences_letters_segment$order_update <- rss_sequences_letters_segment$order
  rss_sequences_count_segment$order_update <- rss_sequences_count_segment$order
  rss_sequences_tag_segment$order_update <- rss_sequences_tag_segment$order
  
  rss_sequences_letters_segment[gene_type != larger_group_id, order_update := order_update + diff]
  rss_sequences_count_segment[gene_type != larger_group_id, order_update := order_update + diff]
  rss_sequences_tag_segment[gene_type != larger_group_id, order_update := order_update + diff]
  
  rss_sequences_tag_segment[,letter_underline:=ifelse(
    underline == "underline",
    sprintf("*%s",letter),
    letter
  )]
  rss_sequences_tag_segment[, letter_underline := gsub("[.][0-9]+$", "", letter_underline)]
  
  p_tag <- ggplot() +
    geom_text_repel(
      data = rss_sequences_tag_segment,
      mapping = aes(
        as.numeric(position), order_update, 
        label = letter_underline, 
        color = color, hjust = -.1
      ), 
      size = 8, 
      bg.r = .1, 
      force = 0, 
      parse = FALSE  # Avoid interpreting * as a math operator
    ) + 
    scale_y_continuous(limits = c(-0.55, max_height_larger_group)) +
    scale_color_manual(values = scale_color_vector) +
    theme_logo() + 
    theme(
      legend.position = 'none', 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0.1, 0, "cm"),
      axis.ticks.length.x=unit(0,'lines'),
      panel.spacing = unit(0, "cm")
    )
  
  
  p_sequences <- ggplot() +
    geom_text(data = rss_sequences_letters_segment,
                    mapping = aes(
                      x = as.numeric(position), 
                      y = order_update, 
                      label = gsub("N", " ", letter), 
                      color = color, 
                      size = ifelse(letter == ".", 16, 8)
                    )) +
    scale_x_continuous(breaks = 1:17, expand = c(0.03, 0)) + 
    scale_y_continuous(limits = c(-0.55, max_height_larger_group)) +
    scale_color_manual(values = scale_color_vector) +
    scale_size_identity() + 
    theme_logo() + 
    theme(
      legend.position = 'none', 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.length.x=unit(0,'lines'),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      panel.spacing = unit(0, "cm")
    )
  

  
  
  p_count <- ggplot() +
    geom_text_repel(data = rss_sequences_count_segment, 
                    mapping = aes(as.numeric(position), order_update, 
                                  label=letter, color=color, hjust = 0,
                                  size = ifelse(letter == ".", 16, 8)),
                    bg.r = .1, 
                    force = 0)+
    scale_color_manual(values=scale_color_vector) +
    scale_y_continuous(limits = c(-0.55, max_height_larger_group)) +
    theme_logo() + 
    scale_size_identity() +
    theme(legend.position = 'none', 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x=unit(0,'lines'),
          plot.margin = margin(0, -0.5, 0, 0, "cm"),
          panel.spacing = unit(0, "cm")
    )
  
  if(nrow(rss_sequences_macaque[gene_type==gene_type_,])>15){
    loc <- max(rss_sequences_letters_segment$order_update)/2
    
  }else{
    loc <- length(unique(rss_sequences_letters_segment$rss_tag))/2
    loc <- max(rss_sequences_letters_segment$order_update)-loc
  }
  
  p_sequences <- p_sequences +
    annotate("text", x = 9, y = loc, label = lab, size = 16)
  
  
  
  
  title <- gene_type_
  if(grepl("D", title)){
    title <- strsplit(title, "_")[[1]]
    title <- paste0(title[2], "' ", title[1])
  }
  
  p_seqlogo <- ggplot() + 
    geom_logo(rss_sequences_data_segment$rss_sequence, 
              method = "probability", seq_type = "dna") +
    theme_logo() + 
    theme(legend.position = 'none', 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = .5, size = 40),
          axis.ticks.length.x=unit(0,'lines')) +
    labs(title = title) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt")) +
    scale_x_continuous(breaks = 1:19,expand = c(0, 0)) 
    
  
  layout <- "
              #A#
              BCD
              "
  
  p <- wrap_elements(
    p_seqlogo + 
      p_tag + 
      p_sequences + 
      p_count + 
      plot_layout(design = layout, widths = c(1/3,1,1/5), heights = c(1/2,1))
  )
  
  rss_by_gene_plot_list[[gene_type_]] <- p
  
}


### add here the diversity.


p <- wrap_elements(wrap_plots(rss_by_gene_plot_list, ncol = 4))

df <- fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_set_with_novel_2025-02-05_with_asc.csv")
df <- df %>%
  dplyr::filter(sample_count_AIRRseq > 0, sure_subject_count>0)

# Summarize the number of unique alleles per gene type
summary_df <- df %>%
  dplyr::group_by(gene_type) %>%
  dplyr::summarize(unique_allele_count = n_distinct(allele))

summary_df <- summary_df %>%
  dplyr::mutate(gene_type = ifelse(gene_type == "IGHD", "IGHD_3", gene_type)) %>%
  dplyr::bind_rows(data.frame(gene_type = "IGHD_5", unique_allele_count = summary_df$unique_allele_count[summary_df$gene_type == "IGHD"]))

summary_df$gene_type <- factor(summary_df$gene_type, levels = c("IGHV","IGHJ","IGLV","IGLJ","IGKV","IGKJ","IGHD_3","IGHD_5"))



hill_diversity <- function(q, proportions) {
  if (q == 1) {
    return(exp(-sum(proportions * log(proportions), na.rm = TRUE)))
  } else {
    return((sum(proportions^q))^(1 / (1 - q)))
  }
}

# Define q values for Hill diversity calculation
q_values <- seq(0, 4, by = 0.1)


# Initialize an empty list to store diversity data for each gene type
diversity_data <- list()

# Loop through each gene type and calculate diversity values
for (gene_type_ in c("IGHV","IGHJ","IGLV","IGLJ","IGKV","IGKJ","IGHD_3","IGHD_5")) {
  # Subset data for the current gene type
  data <- rss_sequences_macaque[gene_type == gene_type_, .(rss_aligned, count)]
  data$proportion <- data$count / sum(data$count)
  
  # Calculate diversity values for each q value
  diversity_values <- sapply(q_values, hill_diversity, proportions = data$proportion)
  
  # Combine q values and diversity values into a data frame
  diversity_df <- data.frame(q = q_values, Diversity = diversity_values, gene_type = gene_type_)
  
  # Append to the diversity data list
  diversity_data[[gene_type_]] <- diversity_df
}

# Combine all gene type diversity data into a single data frame
all_diversity_data <- do.call(rbind, diversity_data)


all_diversity_data <- do.call(rbind, diversity_data) %>%
  left_join(summary_df, by = "gene_type")


all_diversity_data$gene_type <- factor(
  all_diversity_data$gene_type, 
  levels = c("IGHV","IGHJ","IGLV","IGLJ","IGKV","IGKJ","IGHD_3","IGHD_5")
)

first_rss_values <- all_diversity_data %>%
  dplyr::filter(q == 0) %>%
  dplyr::select(gene_type, Diversity) %>%
  dplyr::rename(RSS_count = Diversity)

# Create the initial p_div with a log2 y-axis scale
p_div <- ggplot(all_diversity_data, aes(x = q, y = Diversity)) +
  geom_line() +
  labs(x = "q", y = expression(''^q * D)) +
  theme_minimal() +
  facet_wrap(~ gene_type, nrow = 4, labeller = as_labeller(
    c(
      "IGHV" = "IGHV",
      "IGLV" = "IGLV",
      "IGKV" = "IGKV",
      "IGHD_3" = "3' IGHD",
      "IGHJ" = "IGHJ",
      "IGLJ" = "IGLJ",
      "IGKJ" = "IGKJ",
      "IGHD_5" = "5' IGHD"
    )
  )) +
  theme(
    legend.position = "none",              # Remove legend
    axis.title = element_text(size = 28),  # Increase axis title size
    axis.text = element_text(size = 28),   # Increase axis text size
    strip.text = element_text(size = 28)#,  # Increase facet label text size
    #panel.spacing = unit(5, "lines")
  ) +
  scale_y_continuous(trans = "log2", expand = expansion(mult = c(0.1, 0.2)))

# Add the unique allele counts as a bar plot to the line plot facets
p_div <- p_div + 
  geom_col(data = summary_df, aes(x = min(all_diversity_data$q) - 1, y = unique_allele_count),
           fill = "lightblue", width = 1, inherit.aes = FALSE) +
  geom_text(data = summary_df, aes(x = min(all_diversity_data$q) - 1, y = unique_allele_count, label = unique_allele_count),
            vjust = -0.5, size = 7, inherit.aes = FALSE)  # Text above each bar

# Add arrows and RSS count labels
p_div <- p_div +
  geom_segment(data = first_rss_values %>%
                 left_join(summary_df, by = "gene_type"),
               aes(x = min(all_diversity_data$q) + 1, y = RSS_count + 15,
                   xend = -0.1, yend = RSS_count),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"), # Bold arrow with closed end
               size = 4, inherit.aes = FALSE) +
  geom_text(data = first_rss_values,
            aes(x = min(all_diversity_data$q) + 1, y = RSS_count + 15, label = RSS_count),
            vjust = -1, size = 7, inherit.aes = FALSE)  # Place RSS count label


heptamer_nonamer <- rss_sequences_macaque[,.(count = sum(count)), by=.(gene_type,rss_sequence)]

heptamer_nonamer_wider <- heptamer_nonamer %>% select(gene_type, rss_sequence) %>%
  dplyr::mutate(presence = TRUE) %>% # Add a column to indicate presence
  pivot_wider(names_from = gene_type, values_from = presence, values_fill = FALSE)

data <- heptamer_nonamer_wider
data2 <- heptamer_nonamer

source("scripts/upset_ComplexUpset_v2.R")

p_upset <- upset_adapted_v2(heptamer_nonamer_wider, heptamer_nonamer, 
                         intersect = c("IGHV","IGLV","IGKV","IGHD_3","IGHJ","IGLJ","IGKJ","IGHD_5"),
                         name="Heptamer-Nonamer",
                         themes=upset_modify_themes(
                           list(
                             'intersections_matrix'=theme(text=element_text(size=32), 
                                                          axis.title.x = element_blank(),
                                                          axis.text.x = element_blank(),
                                                          axis.ticks.length.x=unit(0,'lines'),
                                                          panel.grid.minor.y = element_blank()),
                             'overall_sizes'=theme(axis.text.x=element_text(size=32),
                                                   axis.title.x = element_text(size=32))
                           )
                         ),
                         base_annotations=list(
                           'Intersection\nsize'=intersection_size2(
                             counts=TRUE,
                             mapping=aes(fill=rss_group),
                             text_colors = c("black","black"), text = list(size = 8)
                           )+theme(axis.text.y=element_text(size=32),
                                   axis.title.y = element_text(size=32),
                                   legend.text = element_text(size=20),
                                   legend.title = element_text(size=21)) +
                             guides(fill="none") +
                             labs(fill="") + coord_cartesian(clip = "off")
                         ), height_ratio=6)

p_upset_leg <- wrap_elements(p_upset$plot + inset_element(p_upset$leg, -0.5, 1.2, 0, 1,align_to = 'full', on_top = F, clip = T))
# Combine p and p_div, ensuring alignment of panels
layout <- "
AA
BC
"

p_rss <-  wrap_elements(wrap_plots(rss_by_gene_plot_list, ncol = 4) & theme(plot.margin = margin(0, 0, -10, 0, "lines")))

combined_plot <- free(p_rss) + p_div + p_upset_leg + 
  plot_layout(heights = c(2, 1), 
              widths = c(1,2), design = layout) +  # Adjust layout for alignment
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 70))

# Save the final combined plot
ggsave("figures/figure6.pdf", combined_plot, width = 45, height = 45, 
       limitsize = FALSE, units = "in", dpi = 800)
