### script for figure 4

# Function to propagate parent-child relationships and stop when encountering another family
propagate_ancestors <- function(family, p_data, families) {
  
  # Get the child nodes that belong to the current family
  labs <- names(families)[families == family]
  fam_data <- as_tibble(p_data) %>%
    dplyr::filter(label %in% labs)
  
  # Initialize a vector to store the family nodes (child nodes at the beginning)
  all_family_nodes <- fam_data$node
  
  # Track whether new parent nodes were found in the previous step
  new_parents_found <- TRUE
  
  while (new_parents_found) {
    # Find the parent nodes of the current family nodes
    parent_nodes <- p_data %>%
      dplyr::filter(node %in% all_family_nodes) %>%
      dplyr::pull(parent) %>%
      unique()
    
    # Filter for parents where all children belong to the same family
    valid_parents <- p_data %>%
      dplyr::filter(parent %in% parent_nodes) %>%
      dplyr::group_by(parent) %>%
      dplyr::filter(all(node %in% all_family_nodes)) %>%
      dplyr::pull(parent) %>%
      unique()
    
    valid_parents <- valid_parents[! valid_parents %in% all_family_nodes]
    # Check if new valid parents were found
    if (length(valid_parents) > 0) {
      all_family_nodes <- unique(c(all_family_nodes, valid_parents))
    } else {
      new_parents_found <- FALSE
    }
  }
  
  return(all_family_nodes)
}

if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  data.table,
  pbapply,
  ggplot2,
  patchwork,
  ggrepel,
  igraph,
  cluster,
  fpc,
  ggExtra,
  colorspace,
  pracma,
  ggtree,
  gridtext,
  ggforce,
  dplyr
)

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

load('data/figure4/figure4_data.rda')

MUSA_set <- fread("data/MUSA_data/MUSA_set_all_2025-02-05_with_asc.csv")


pal_family <- c(
  hcl.colors(12, palette = "Zissou1"),
  hcl.colors(12, palette = "Roma"),
  hcl.colors(12, palette = "Spectral"),
  hcl.colors(12, palette = "Cividis"),
  hcl.colors(12, palette = "Temps"),
  hcl.colors(12, palette = "Earth")
)

macaque_segments <- c("macaque_IGHV", "macaque_IGKV", "macaque_IGLV", "macaque_IGHJ", "macaque_IGKJ", "macaque_IGLJ", "macaque_IGHD")

max_fam <- unique(sort(unlist(lapply( macaque_segments, function(segment){
  chain <- gsub("macaque_","",segment)
  as.numeric(gsub(chain,"",unique(alakazam::getFamily(clusters_tables[[segment]]$new_name, strip_d = F, omit_nl = F))))
}))))

pal_families <- setNames(pal_family[1:length(max_fam)], as.character(max_fam))

pal_families <- c(pal_families, setNames("gray80",0))

tree_metadata_all_segments <- list()
segments_wheels <- list()
for(segment in macaque_segments){
  iglabel_new_names <- setNames(clusters_tables[[segment]]$new_tag, clusters_tables[[segment]]$original_allele)
  
  ref_dist_segment <- ref_distance[[segment]]$matrix
  colnames(ref_dist_segment) <- iglabel_new_names[colnames(ref_dist_segment)]
  rownames(ref_dist_segment) <- iglabel_new_names[rownames(ref_dist_segment)]
  
  metadata_segment <- setDT(clusters_tables[[segment]])
  
  if(grepl("IGHD",segment)){
    ## for the figure we will filter out the non-sure IGHD novel alleles. 
    alleles <- MUSA_set[gene_type=="IGHD" & novel==TRUE & sample_count_AIRRseq==0,new_tag]
    metadata_segment <- metadata_segment[!new_tag %in% alleles]
    ref_dist_segment <- ref_dist_segment[!rownames(ref_dist_segment) %in% alleles,!colnames(ref_dist_segment) %in% alleles]
  }
  
  
  ref_dist_segment <- as.dist(ref_dist_segment)
  
  
  
  hc <- hclust(ref_dist_segment, method = "complete")
  dend <- as.dendrogram(hc, hang = -1)
  
  tree <- ape::as.phylo(hc)
  
  asc_tips <- tapply(metadata_segment$new_tag, metadata_segment$new_cluster, function(x) x)
  
  collapse_clade <- function(tree, node) {
    ape::drop.tip(tree, tip = node[-1])
  }
  
  collapsed_tree <- Reduce(function(tree, node) collapse_clade(tree, node), asc_tips, init = tree)
  
  hc_collapsed_tree <- ape::as.hclust.phylo(collapsed_tree)
  
  p_asc <- ggtree(hc_collapsed_tree, layout='circular', size=2) 
  
  p_asc_dendo <- ggtree(hc_collapsed_tree, layout='dendrogram', size=2)
  
  p_data <- p_asc$data
  p_data_dendo <- p_asc_dendo$data
  
  p_data$branch_dendo <- p_data_dendo$branch
  
  labels <- p_data$label
  labels <- labels[!is.na(labels)]
  
  families <- alakazam::getFamily(labels, strip_d = F, omit_nl = F)
  families <- setNames(families, labels)
  
  tree_metadata <- data.frame(node = integer(),
                              label = character(),
                              family = character(),
                              asc = character(),
                              size = numeric(),
                              BG505 = logical()
  )
  
  # Loop through each family and propagate ancestors
  tree_metadata <- rbindlist(lapply(unique(families), function(fam){
    # Propagate ancestors for each family
    fam_nodes <- propagate_ancestors(family = fam, p_data = as_tibble(p_data), families = families)
    # Add the results to the tree metadata
    lab <- p_data$label[p_data$node %in% fam_nodes]
    present <- sapply(lab, function(l){
      if(!is.na(l)){
        l <- strsplit(l,"[*]")[[1]][1]
        any(metadata_segment[grepl(l,new_tag), present])
      }else{
        FALSE
      }
    }, USE.NAMES = F)
    
    data.frame(
      node = fam_nodes,
      label = lab,
      present = present,
      family = fam,
      asc = "",
      size = NA
    )
  }))
  
  tree_metadata_dendo <- rbindlist(lapply(unique(families), function(fam){
    # Propagate ancestors for each family
    fam_nodes <- propagate_ancestors(family = fam, p_data = as_tibble(p_data_dendo), families = families)
    # Add the results to the tree metadata
    lab <- p_data_dendo$label[p_data$node %in% fam_nodes]
    present <- sapply(lab, function(l){
      if(!is.na(l)){
        l <- strsplit(l,"[*]")[[1]][1]
        any(metadata_segment[grepl(l,new_tag), present])
      }else{
        FALSE
      }
    }, USE.NAMES = F)
    
    data.frame(
      node = fam_nodes,
      label = lab,
      present = present,
      family = fam,
      asc = "",
      size = NA
    )
  }))
  
  chain <- substr(gsub("macaque_","",segment),1,3)
  
  tree_metadata$asc[!is.na(tree_metadata$label)] <- alakazam::getGene(tree_metadata$label[!is.na(tree_metadata$label)], strip_d = F, omit_nl = F)
  tree_metadata$asc_label[!is.na(tree_metadata$label)] <- gsub(chain,"",tree_metadata$asc[!is.na(tree_metadata$label)])
  
  tree_metadata_dendo$asc[!is.na(tree_metadata_dendo$label)] <- alakazam::getGene(tree_metadata_dendo$label[!is.na(tree_metadata_dendo$label)], strip_d = F, omit_nl = F)
  tree_metadata_dendo$asc_label[!is.na(tree_metadata_dendo$label)] <-  gsub(chain,"",tree_metadata_dendo$asc[!is.na(tree_metadata_dendo$label)])
  
  ## create a metadata of the number of alleles per per node
  for(a in labels){
    g <- metadata_segment[new_tag == a, new_cluster]
    size <- nrow(metadata_segment[new_cluster == g])
    tree_metadata$size[tree_metadata$label == a] <- size
  }
  
  pal <- setNames(hcl.colors(length(unique(families)), palette = "Zissou1"), unique(families))
  
  cut_height <- max(nchar(gsub("[.]","",ref[[segment]]))) * 0.25
  
  x = abs(p_data_dendo$branch)
  p_data_dendo$branch_dendo_abs <- abs(p_data_dendo$branch)
  p_data_dendo$branch_dendo_dist <- abs(p_data_dendo$branch_dendo_abs- cut_height)
  cut_height_dendo <- which(abs(x - cut_height) == min(abs(x - cut_height)))
  cut_height_dendo <- cut_height_dendo[which.min(p_data_dendo$branch[cut_height_dendo])]
  cut_height_dendo_diff <- -p_data_dendo$branch_dendo_dist[cut_height_dendo]
  
  # get the radius for the 75% threshold
  radius_cutoff <- p_data$branch[cut_height_dendo] + cut_height_dendo_diff
  
  tree_metadata <- tree_metadata %>%
    mutate(size_bin = cut(size, 
                          breaks = 2^(0:ceiling(log2(max(size, na.rm = TRUE)))), 
                          include.lowest = TRUE))
  
  
  bin_labels <- levels(tree_metadata$size_bin)
  
  sizes <- seq(5, 15, by=1.5)[1:length(bin_labels)]
  
  label_size <- 7
  
  
  tree_metadata$family_tag <- as.numeric(gsub(gsub("macaque_","",segment),"",tree_metadata$family))
  
  offset <- floor(log2(length(labels)))
  offset <- offset*ifelse(segment=="macaque_IGHV", 0.5, ifelse(grepl("V",segment), 0.85, 0.85))
  
  tree_metadata_all_segments[[segment]] <- tree_metadata
  segments_wheels[[segment]] <- p_asc %<+% tree_metadata + aes(color = as.factor(family_tag)) + 
    scale_color_manual(values = pal_families, breaks = ~ .x[!is.na(.x)]) +
    geom_tippoint(aes(size=as.factor(size_bin), fill=present), alpha=0.8, color="gray90", shape=21)+#color="gray40") +
    scale_fill_manual(values=c("TRUE" = "tomato2","FALSE" = "gray40"), breaks = c("TRUE" = "Present")) +
    ggnewscale::new_scale_colour() + 
    geom_tiplab(aes(label=asc_label), show.legend=FALSE, offset=offset, size=label_size, color = "black") +
    labs(size="# Alleles") + 
    guides(size = guide_legend(override.aes = list(fill = "gray40"))) +
    geom_vline(xintercept = radius_cutoff, colour = "firebrick", linetype = "dashed") +
    scale_size_manual(values = sizes,
                      labels = bin_labels) +
    theme_void() + 
    theme(legend.position = "none"
    )
  
  if(any(grepl("64",bin_labels))){
    
    allele_legend <- p_asc %<+% tree_metadata  +
      geom_tippoint(aes(size=as.factor(size_bin)), alpha=0.7, color="gray40") +
      labs(size="# Alleles") +
      guides(size = guide_legend(override.aes = list(fill = "gray40"), theme = theme(
        legend.text = element_text(size=35), legend.background = element_blank(),
        legend.title = element_text(size=35)
      ))) +
      scale_size_manual(values = sizes,
                        labels = bin_labels)
    
    legend <- ggfun::get_legend(allele_legend)
    
  }
  
  
}

tree_metadata_all_segments_df <- rbindlist(tree_metadata_all_segments, idcol = "gene_type")
tree_metadata_all_segments_df[,gene_type:=gsub("macaque_","",gene_type)]
# remove rows that the label is na
tree_metadata_all_segments_df <- tree_metadata_all_segments_df[!is.na(label)]

levels <-  c("IGHV","IGKV","IGLV","IGHD","IGHJ","IGKJ","IGLJ")
tree_metadata_all_segments_df$gene_type <- factor(tree_metadata_all_segments_df$gene_type, levels=levels)

bins <- c(0, 1, 2, 4, 8, 16, 32, 64, Inf)
labels <- c("(0,1]", "(1,2]", "(2,4]", "(4,8]", "(8,16]", "(16,32]", "(32,64]", "(64,128]")
tree_metadata_all_segments_df$bin_size <- cut(tree_metadata_all_segments_df$size, breaks = bins, labels = labels, right = TRUE)
tree_metadata_all_segments_df$present_lab <- ifelse(tree_metadata_all_segments_df$present, "In AIRR-seq", "Not in AIRR-seq")

colors_seen <- setNames(c("plum4","sandybrown"),c( "In AIRR-seq", "Not in AIRR-seq"))

vs <- ggplot(tree_metadata_all_segments_df[grepl("V",gene_type)], 
             aes(x = bin_size, fill = present_lab)) + 
  geom_bar(alpha = 1) + 
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  ) +
  theme_minimal(base_size = 34) + 
  facet_wrap(.~gene_type, nrow = 1, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'top') + labs(x = "# Alleles", y = "# ASCs") +
  guides(fill = "none")

jd <- ggplot(tree_metadata_all_segments_df[!grepl("V",gene_type)], 
             aes(x = bin_size, fill = present_lab)) + 
  geom_bar(alpha = 1) + 
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  ) +
  theme_minimal(base_size = 34) + facet_wrap(.~gene_type, nrow = 1, scales = "free") + guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x = "# Alleles", y = "# ASCs")

p_allele_to_asc <- wrap_elements(vs / jd)

clusters_tables_df <- rbindlist(clusters_tables)

l <- setNames(clusters_tables_df$new_name, clusters_tables_df$new_tag)

tree_metadata_all_segments_df_family <- tree_metadata_all_segments_df[,new_allele:=l[label]][
  ,family:=alakazam::getFamily(new_allele)][,.(size = sum(size)), by = .(gene_type, family, present)]
bins <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
labels <- c("(0,1]", "(1,2]", "(2,4]", "(4,8]", "(8,16]", "(16,32]", "(32,64]", "(64,128]", "(128,256]", "(256,512]", "(512,1024]")
tree_metadata_all_segments_df_family$bin_size <- cut(tree_metadata_all_segments_df_family$size, breaks = bins, labels = labels, right = TRUE)
tree_metadata_all_segments_df_family$present_lab <- ifelse(tree_metadata_all_segments_df_family$present, "In AIRR-seq", "Not in AIRR-seq")

vs <- ggplot(tree_metadata_all_segments_df_family[grepl("V",gene_type)], 
             aes(x = bin_size, fill = present_lab)) + 
  geom_bar(alpha = 1) + 
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  ) +
  theme_minimal(base_size = 34) + 
  facet_wrap(.~gene_type, nrow = 1, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'top') + labs(x = "# Alleles", y = "# Families") + guides(fill = "none")

jd <- ggplot(tree_metadata_all_segments_df_family[!grepl("V",gene_type)], 
             aes(x = bin_size, fill = present_lab)) + 
  geom_bar(alpha = 1) + 
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  ) +
  theme_minimal(base_size = 34) + facet_wrap(.~gene_type, nrow = 1, scales = "free") + guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x = "# Alleles", y = "# Families")

p_allele_to_family <- wrap_elements(vs / jd)


## to bar plots per each gene type number of families and number of asc
clusters_tables_df[,gene_type:=substr(new_name,1,4)]
clusters_tables_df[,asc:=alakazam::getGene(new_name, strip_d = F, omit_nl = F)]
clusters_tables_df[,family:=alakazam::getFamily(new_name, strip_d = F, omit_nl = F)]

## filter the d alleles
d_alleles <- MUSA_set[gene_type=="IGHD" & novel==TRUE & sample_count_AIRRseq==0,new_tag]
asc_summary <- clusters_tables_df[!new_tag %in% d_alleles,]
asc_summary <- asc_summary[,.(alleles = .N, present = any(present)), by = .(gene_type, asc, family)]
asc_summary <- asc_summary[,.(asc = .N, mean_alleles = mean(alleles), present_asc = sum(present), present = any(present)), by = .(gene_type, family)]
asc_summary_all <- asc_summary[,.(family = .N, asc = sum(asc), present_asc = sum(present_asc), present_family = sum(present)), by = .(gene_type)]
asc_summary_all[, not_present_asc := asc - present_asc] # ASCs not present
asc_summary_all[, not_present_family := family - present_family] # Families not present
# Reshape the data to long format for stacked barplot
barplot_table <- data.table(
  gene_type = rep(asc_summary_all$gene_type, 4),
  type = rep(c("Families", "Families", "ASCs", "ASCs"), each = nrow(asc_summary_all)),
  presence = rep(c("In AIRR-seq","Not in AIRR-seq", "In AIRR-seq","Not in AIRR-seq"), each =  nrow(asc_summary_all)),
  count = c(
    asc_summary_all$present_family, asc_summary_all$family - asc_summary_all$present_family,
    asc_summary_all$present_asc, asc_summary_all$asc - asc_summary_all$present_asc
  )
)

levels <-  c("IGHV","IGKV","IGLV","IGHD","IGHJ","IGKJ","IGLJ")
barplot_table$gene_type <- factor(barplot_table$gene_type, levels=levels)

vs <- ggplot(barplot_table[grepl("V",gene_type)], aes(x = type, y = count, fill = presence)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(.~gene_type, scales = "free") +
  labs(x = "",y = "Count",fill = "") +
  theme_minimal(base_size = 34) +
  guides(fill = "none") +
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  )

jd <- ggplot(barplot_table[!grepl("V",gene_type)], aes(x = type, y = count, fill = presence)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(.~gene_type, scales = "free") +
  labs(x = "",y = "Count",fill = "") +
  theme_minimal(base_size = 34) +
  scale_fill_manual(
    values = colors_seen,
    drop = FALSE,
    name = ""
  ) + theme(legend.position = 'top', legend.text = element_text(size=35), legend.background = element_blank(),
            legend.title = element_text(size=35),
            legend.key = element_blank())

legend_seen <- ggfun::get_legend(jd)

for (i in seq_along(legend_seen$grobs)) {
  if (inherits(legend_seen$grobs[[i]], "rect")) {  # Check if the grob is a rectangle (background element)
    legend_seen$grobs[[i]]$gp <- gpar(fill = "transparent", col = NA)  # Make background transparent
  }
}

jd <- jd + guides(fill = "none")

p_asc_family <- wrap_elements(vs / jd)


layout2 <- "
ABD
ABE
AHH
ACF
ACG
"

patch <- 
  (segments_wheels[["macaque_IGHV"]] + 
     scale_x_continuous(limits = c(0, NA), # set limit the same as in the y axis https://github.com/YuLab-SMU/ggtree/blob/c17773c973d6c4036ee3af40a3957fb74d8ee9ff/R/method-ggplot-add.R#L110
                        expand = ggplot2::expansion(add = c(0, -40)) # use the add parameter in the ggplot2::ggplot2::expansion function instead of mult. This will expand the x-scale without altering the center of the dendrogram.
                        # the addition value was set based on optimization to the pdf, it is not a general rule.
     ) + ggtitle("IGHV") + theme(plot.title = element_text(hjust = 0.5, vjust = 6, size = 45)) +
     inset_element(legend, 0.15, 0.98, 0, 1, align_to = "panel", on_top = F, clip = T)) +  
  (segments_wheels[["macaque_IGKV"]] + 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, -8))) + 
     ggtitle("IGKV") + theme(
       plot.title = element_text(hjust = 0.5, vjust = 2.5, size = 45),
       plot.background = element_rect(color = "blue", linewidth = 2))) +
  (segments_wheels[["macaque_IGLV"]] + 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, -8))) + 
     ggtitle("IGLV") + theme(
       plot.title = element_text(hjust = 0.5, vjust = 2, size = 45),
       plot.background = element_rect(color = "red", linewidth = 2))) +
  (rotate_tree(segments_wheels[["macaque_IGHD"]], 75)+ 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, 15)))+
     ggtitle("IGHD") + theme(
       plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 45),
       plot.background = element_rect(color = "brown", linewidth = 2))) +
  (rotate_tree(segments_wheels[["macaque_IGHJ"]], 75) + 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, 8)))+
     ggtitle("IGHJ") + theme(
       plot.title = element_text(hjust = 0.5, vjust = -1, size = 45),
       plot.background = element_rect(color = "orange", linewidth = 2))) + 
  (rotate_tree(segments_wheels[["macaque_IGKJ"]], 75) + 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, 8)))+
     ggtitle("IGKJ") + theme(
       plot.title = element_text(hjust = 0.5, vjust = -1, size = 45),
       plot.background = element_rect(color = "green", linewidth = 2))) +  
  (rotate_tree(segments_wheels[["macaque_IGLJ"]], 75) + 
     scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(add = c(0, 8)))+
     ggtitle("IGLJ") + theme(
       plot.title = element_text(hjust = 0.5, vjust = -1, size = 45),
       plot.background = element_rect(color = "purple", linewidth = 2))) + 
  plot_spacer() + 
  plot_layout(design = layout2, 
              widths = c(4,2,1), # column widths, main = 4, sub1 = 2, sub2 = 1
              heights = c(6,6,0.2,6,6) # all rows have the same heights
  ) & theme(plot.margin = margin(0, 0, 0.5, 0, "lines"), # creating a small space between the columns, r=2
            plot.background = element_blank()
  )

patch_wrap <- wrap_elements(patch)

layout_final <- "
DDD
ABC
"
figure_3 <-  
  p_allele_to_asc + 
  wrap_elements(p_allele_to_family + inset_element(legend_seen, 1, 1.1, 0, 1, 
                                                   align_to = "panel", on_top = T, clip = F, ignore_tag = T))+ 
  p_asc_family + 
  (patch_wrap) + 
  plot_layout(heights = c(2,1), design = layout_final) +
  plot_annotation(tag_levels = list(c("B","C","D","A"))) & 
  theme(plot.tag = element_text(size = 80, hjust = -1, vjust = 0),
        plot.margin = margin(0, 0.1, 0, 0.1, "lines"))


ggsave(filename = "figures/figure4.pdf", 
       figure_3, width = 50, height = 45, units = "in", limitsize = F)
