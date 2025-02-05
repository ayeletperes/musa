
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
library(ggplot2)
library(VennDiagram)
library(grid)
library(patchwork)
library(cowplot)
library(Biostrings)


IGH_D <- readDNAStringSet("SeedSet/2025-01-03/IGH/D.fasta")
IGH_V <- readDNAStringSet("SeedSet/2025-01-03/IGH/V_gapped.fasta")
IGH_J <- readDNAStringSet("SeedSet/2025-01-03/IGH/J.fasta")
IGL_V <- readDNAStringSet("SeedSet/2025-01-03/IGL/V_gapped.fasta")
IGL_J <- readDNAStringSet("SeedSet/2025-01-03/IGL/J.fasta")
IGK_V <- readDNAStringSet("SeedSet/2025-01-03/IGK/V_gapped.fasta")
IGK_J <- readDNAStringSet("SeedSet/2025-01-03/IGK/J.fasta")

# Combine all alleles from FASTA files into a named list
reference_alleles <- list(
  IGH_D = names(IGH_D),
  IGH_V = names(IGH_V),
  IGH_J = names(IGH_J),
  IGL_V = names(IGL_V),
  IGL_J = names(IGL_J),
  IGK_V = names(IGK_V),
  IGK_J = names(IGK_J)
)

all_reference_alleles <- unlist(reference_alleles)

selected_columns <- c(
  "asc", "asc_gene","subject"
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

# Combine data frames
digger_df <- rbind(digger_df_IGH, digger_df_IGL, digger_df_IGK)


digger_df <-  digger_df %>%
 mutate(chain = str_sub(asc, 1,3),
         call = str_sub(asc, 4,4),
         gene_type = str_sub(asc, 1,4))

digger_df <- unique(digger_df)

actual_df <- read.csv("data/rhesus_macaque_data/repertoire_genotype_all_04_02_rep_data.csv")



actual_df <-  actual_df %>%
 mutate(chain = str_sub(allele, 1,3),
         call = str_sub(allele, 4,4),
         gene_type = str_sub(allele, 1,4))
  


personal_color <- "#00519F"
gg_color <-  "grey70"

allele_sample_count <- actual_df %>%
  dplyr::mutate(chain = str_sub(allele, 1, 3),
         call = str_sub(allele, 4, 4),
         gene_type = str_sub(allele, 1, 4)) %>%
  dplyr::group_by(allele, chain, call, gene_type) %>%
  dplyr::summarise(sample_count = n_distinct(sample), .groups = "drop")


cdf_data <- allele_sample_count %>%
  dplyr::group_by(chain, call, gene_type, sample_count) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%  # Count occurrences of each sample_count
  dplyr::arrange(sample_count, .by_group = TRUE) %>%  # Ensure ascending order
  dplyr::group_by(chain, call, gene_type) %>%
  dplyr::mutate(cumulative_alleles = cumsum(count)) %>%  # Calculate cumulative sum
  dplyr::ungroup()


cdf_data <- cdf_data %>%
  dplyr::mutate(source = "ERP")  # First dataset



allele_sample_count_gg <- digger_df %>%
 dplyr::mutate(chain = str_sub(asc, 1,3),
         call = str_sub(asc, 4,4),
         gene_type = str_sub(asc, 1,4)) %>%
  dplyr::group_by(asc, chain, call,gene_type) %>%
  dplyr::summarise(sample_count = n_distinct(subject))

cdf_data_gg <- allele_sample_count_gg %>%
  dplyr::group_by(chain, call, gene_type, sample_count) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%  # Count occurrences of each sample_count
  dplyr::arrange(sample_count, .by_group = TRUE) %>%  # Ensure ascending order
  dplyr::group_by(chain, call, gene_type) %>%
  dplyr::mutate(cumulative_alleles = cumsum(count)) %>%  # Calculate cumulative sum
  dplyr::ungroup()

cdf_data_gg <- cdf_data_gg %>%
  dplyr::mutate(source = "GG")  # Sec dataset


# Combine the datasets and sep
combined_data <- bind_rows(cdf_data, cdf_data_gg)
erp_data <- combined_data %>% filter(source == "ERP")
gg_data <- combined_data %>% filter(source == "GG")


create_cdf_chain <- function(chain_name) {
    x <- erp_data[erp_data$gene_type==chain_name,]
    t <- gg_data[gg_data$gene_type==chain_name,]

    scaleFactor <- max(x$cumulative_alleles) / max(t$cumulative_alleles)

  theme_custom <- if (chain_name == "IGKJ") {
    theme(
      axis.title.y = element_text(color = personal_color, size = 34),
      axis.text.y = element_text(color = personal_color, size = 24),
      axis.title.y.right = element_text(color = gg_color, size = 34),
      axis.text.y.right = element_text(color = gg_color, size = 24),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 10, b = 10, l = 10, r = 10)
    )
  } else if (chain_name == "IGLV") {
    theme(
      axis.text.x = element_text(size = 24),
      axis.title.x = element_text(size = 34),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = personal_color, size = 24),
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_text(color = gg_color, size = 24),
      plot.margin = margin(t = 10, b = 10, l = 10, r = 10)
    )
  } else {
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = personal_color, size = 24),
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_text(color = gg_color, size = 24),
      plot.margin = margin(t = 10, b = 10, l = 10, r = 10)
    )
  }


  get_y_at_x <- function(data, x_value, y_column) {
    if (x_value %in% data$sample_count) {
      # Exact match exists
      return(data[[y_column]][data$sample_count == x_value])
    } else {
      # Interpolate between points if exact match does not exist
      data <- data[order(data$sample_count), ]  # Ensure data is sorted by x
      lower <- max(data$sample_count[data$sample_count < x_value], na.rm = TRUE)
      upper <- min(data$sample_count[data$sample_count > x_value], na.rm = TRUE)
      lower_y <- data[[y_column]][data$sample_count == lower]
      upper_y <- data[[y_column]][data$sample_count == upper]
      return(lower_y + (upper_y - lower_y) * ((x_value - lower) / (upper - lower)))
    }
  }

  # Calculate y-values at x = 10 for both datasets
  erp_yintercept <- get_y_at_x(x, 10, "cumulative_alleles")
  gg_yintercept <- get_y_at_x(t, 10, "cumulative_alleles") * scaleFactor

  # Plot with exact match for y-values at x = 10
  plot <- ggplot() +
      # Primary y-axis for ERP
      geom_line(data = x, aes(x = sample_count, y = cumulative_alleles, color = "ERP"), size = 1.5) +
      scale_y_continuous(
        name = "Expressed - Number of Alleles",
        sec.axis = sec_axis(~./scaleFactor, name = "GG - Number of Alleles")  # Secondary y-axis for GG
      ) +
      # Secondary y-axis for GG
      geom_line(data = t, aes(x = sample_count, y = cumulative_alleles * scaleFactor, color = "GG"), size = 1.5) +
      # Add vertical dashed line at 10 subjects
      geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
      # Add horizontal dashed lines for ERP and GG
      #geom_hline(yintercept = erp_yintercept, linetype = "dashed", color = "blue") +
      #geom_hline(yintercept = gg_yintercept, linetype = "dashed", color = "gray") +
        geom_segment(aes(x = 0, xend = 10, y = erp_yintercept, yend = erp_yintercept), linetype = "dashed", color = personal_color) +
      geom_segment(aes(x = 10, xend = 106, y = gg_yintercept, yend = gg_yintercept), linetype = "dashed", color = gg_color) +
      scale_color_manual(
        values = c("ERP" = personal_color, "GG" = gg_color),  # Assign colors to sources
        guide = "none"  # Remove legend
      ) +
      labs(x = "Number of subjects") +
      theme_minimal() +
      theme_custom

      return(plot)
}

# Create the Venn diagrams for each chain
venn_IGHD <- create_cdf_chain("IGHD")
venn_IGHJ <- create_cdf_chain("IGHJ")
venn_IGHV <- create_cdf_chain("IGHV")
venn_IGLJ <- create_cdf_chain("IGLJ")
venn_IGLV <- create_cdf_chain("IGLV")
venn_IGKJ <- create_cdf_chain("IGKJ")
venn_IGKV <- create_cdf_chain("IGKV")


# Combine all Venn diagrams without duplicating the legend and ensuring consistent heights
plot_cdf <- (venn_IGHD / venn_IGHJ / venn_IGHV / venn_IGKJ / venn_IGKV / venn_IGLJ / venn_IGLV) +
  plot_layout(ncol = 1, heights = c(1, 1, 1,1,1,1,1)) +  # Ensuring equal height for each plot
  theme(
    plot.margin = margin(t = 10, b = -5, l =0, r = 0),  # Adjust margins for better alignment
    panel.spacing = unit(1, "lines")  # Ensure consistent panel spacing
  )


combined_df <- full_join(
  actual_df %>%
    select(allele, sample, in_genomic, chain) %>%
    dplyr::mutate(in_repertior = TRUE),  # Mark rows from actual_df
  digger_df %>%
    dplyr::rename(allele = asc, sample = subject) %>%
    select(allele, sample, chain),
  by = c("allele", "sample", "chain")
)

combined_df <- combined_df %>%
  dplyr::mutate(
    in_repertior = ifelse(is.na(in_repertior), FALSE, in_repertior),  # If NA, mark as FALSE
    in_genomic = ifelse(is.na(in_genomic), TRUE, in_genomic)         # If NA, mark as FALSE
  )

combined_df <- unique(combined_df)

allele_count_per_sample_repertior_not_divid <- combined_df %>%
  dplyr::mutate(call = str_sub(allele, 4,4),
         gene_type = str_sub(allele, 1,4)) %>%
  dplyr::group_by(sample, chain, call,gene_type) %>%
  dplyr::summarise(allele_count = n_distinct(allele))

# Count the number of unique alleles per sample and chain, grouped by call
allele_count_per_sample_repertior <- combined_df %>%
  dplyr::mutate(call = str_sub(allele, 4,4),
         gene_type = str_sub(allele, 1,4)) %>%
  dplyr::group_by(sample, chain, in_repertior, call,gene_type) %>%
  dplyr::summarise(allele_count = n_distinct(allele))



# Step 1: Summarize the allele counts for in_repertior == TRUE for each sample and call
allele_count_repertior_true <- allele_count_per_sample_repertior %>%
  filter(in_repertior == TRUE) %>%
  dplyr::group_by(sample, call, chain, gene_type) %>%
  dplyr::summarise(allele_count_true = sum(allele_count), .groups = 'drop')  # One row per sample and call


plot1 <- ggplot(allele_count_per_sample_repertior_not_divid, aes(x = as.numeric(factor(sample)), y = allele_count)) +
  geom_bar(stat = "identity", fill = gg_color, color = "white") +  # Set the width of the bars to span the entire space
  # Add points for in_repertior == TRUE
  geom_point(
    data = allele_count_repertior_true,
    aes(x = as.numeric(factor(sample)),  # x-value for the point
        y = allele_count_true,           # y-value for the point
        color = "personal_color"),       # Use the personal color
    size = 3,                            # Set the size of the points
    shape = 19                           # Solid circle shape
  ) +
  labs(x = "Subject", y = "Number of Alleles") + # Facet by call within the IGH chain
  theme_minimal() +
  facet_wrap(~gene_type, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_color_manual(values = c("personal_color" = personal_color)) + 
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust =0.5),
    axis.text = element_text(size = 24),          # Increase y-axis text size
    axis.title = element_text(size = 34, margin = margin(r = 40)),  
    strip.text = element_blank(), 
    legend.position = "none",
    plot.margin = margin(t = 10, b = 10, l = 0, r = -5),
    plot.title = element_text(size = 24),
    panel.spacing = unit(2, "lines")
  ) + scale_x_continuous(expand = c(0,0))

plot_strip_only <- ggplot(allele_count_per_sample_repertior, aes(x = as.numeric(factor(sample)), y = allele_count)) +
  geom_blank() +  # No data plotted, just use blank elements
  facet_wrap(
    ~gene_type, 
    scales = "free_y", 
    ncol = 1, 
    strip.position = "left" # Move strip to the left
  ) +
  labs(x = NULL, y = NULL) + # Remove axis labels
  theme_minimal() +
  theme(
    axis.title = element_blank(), # Remove axis titles
    axis.text = element_blank(), # Remove axis text
    axis.ticks = element_blank(), # Remove axis ticks
    panel.grid = element_blank(), # Remove grid lines
    strip.text.y.left = element_text(size = 34, angle = 90, hjust = 0.5, vjust = 0.5), # Keep strip labels
    strip.placement = "outside", # Place strips outside plot area
    plot.margin = margin(t = 10, b = 10, l = 10, r = 10) # Add left margin for alignment
  )


MUSA <- read.csv("MUSA/MUSA_with_rss_and_leader_and_asc_2025-02-05.csv")

#Create a dataframe for alleles in the reference but not in in_repertior_df

reference_df <- data.frame(
  allele = MUSA[MUSA$novel == "FALSE","allele"],
  in_repertior = "FALSE",
  in_reference = "TRUE",
  in_genomic = "FALSE",
  stringsAsFactors = FALSE
)

# Add a new column to in_repertior_df to indicate if the allele is in the references
combined_df$in_reference <- combined_df$allele %in% as.character(all_reference_alleles)

df <- unique(combined_df[,c("allele", "in_repertior", "in_reference", "in_genomic")])

# Combine in_repertior_df and reference_df
combined_with_references <- rbind(
  df, 
  reference_df[!reference_df$allele %in% combined_df$allele, ] # Only add alleles not already in in_repertior_df
)

combined_with_references <- combined_with_references %>%
  dplyr::mutate(chain = str_sub(allele, 1,3),
         gene_type = str_sub(allele, 1,4))

combined_with_references$in_genomic <- ifelse(combined_with_references$in_genomic == "True", TRUE, combined_with_references$in_genomic)



# Step 2: Function to create a Venn diagram with consistent "In Reference" on the left and "In Repertior" on the right

venn_to_ggplot <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
  grob <- grid.grab()
  p <- ggplot() + 
    annotation_custom(grob) + 
    theme_void()  # No axes, titles, etc.
  return(p)
}

create_venn_for_chain <- function(chain_name) {

  chain_data <- combined_with_references %>% filter(gene_type == chain_name)
  
  # Remove NAs from in_repertior and in_reference columns
  chain_data <- chain_data %>% 
    filter(!is.na(in_repertior) & !is.na(in_reference))
  
  # Prepare data for Venn diagram
  venn_data <- list(
    "In Reference" = unique(chain_data$allele[chain_data$in_reference == "TRUE"]),
    "In Repertior" = unique(chain_data$allele[chain_data$in_repertior == "TRUE"]),
    "In Genomic" = unique(chain_data$allele[chain_data$in_genomic == "TRUE"])
  )
  
  # Check for empty sets
  if (length(venn_data$`In Reference`) == 0) {
    venn_data$`In Reference` <- "empty"
  }
  if (length(venn_data$`In Repertior`) == 0) {
    venn_data$`In Repertior` <- "empty"
  }

  rotation.degree = 180
  rev = F
  if(chain_name == "IGKJ"){
    rotation.degree = 0
  }
  if(chain_name %in% c("IGLJ","IGHJ")){
    rev = T
    rotation.degree = 360
  }
  venn_ggplot <- venn_to_ggplot(venn_data, scaled=T, rotation.degree = rotation.degree,
                            col = c("#2ECC71" , personal_color, gg_color), reverse = rev,
                            cat.cex = 0, cex = 2.5, lty = c(1,1,2), lwd = 2) 
  return(venn_ggplot)
}

# Create the Venn diagrams for each chain
venn_IGHD <- create_venn_for_chain("IGHD")
venn_IGHJ <- create_venn_for_chain("IGHJ")
venn_IGHV <- create_venn_for_chain("IGHV")
venn_IGLJ <- create_venn_for_chain("IGLJ")
venn_IGLV <- create_venn_for_chain("IGLV")
venn_IGKJ <- create_venn_for_chain("IGKJ")
venn_IGKV <- create_venn_for_chain("IGKV")


# Combine all Venn diagrams without duplicating the legend and ensuring consistent heights
plot3 <- (venn_IGHD / venn_IGHJ / venn_IGHV / venn_IGKJ / venn_IGKV / venn_IGLJ / venn_IGLV) +
  plot_layout(ncol = 1, heights = c(1, 1, 1,1,1,1,1)) +  # Ensuring equal height for each plot
  theme(
    plot.margin = margin(t = 10, b = -5, l =0, r = 0),  # Adjust margins for better alignment
    panel.spacing = unit(1, "lines")  # Ensure consistent panel spacing
  )


legend_plot <- ggplot(data = data.frame(x = c(1,2,3), y = c(1,2,3), 
                                        label = c("In Repertior", "In Reference", "In Genomic"))) +
  geom_point(aes(x = x, y = y, color = label), size = 4) +  # Invisible points
  scale_color_manual(
    name = "",
    values = c("In Repertior" = personal_color , "In Reference" = "#2ECC71", "In Genomic" = gg_color),
    labels = c("In Repertior" = "Expressed in cohort" , "In Reference" = "Baseline Reference" , "In Genomic" = "GG"), 
    guide = guide_legend(override.aes = list(size = 6), ncol = 3)
  ) +
  theme_void() +
  theme(
    #legend.position = "left",
    #legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.key = element_blank()
  )

legend_plot2 <- ggfun::get_legend(legend_plot)


# Count the number of unique alleles per sample and chain, grouped by call
allele_count_sample <- combined_df %>%
  dplyr::mutate(call = str_sub(allele, 4,4),
         gene_type = str_sub(allele, 1,4)) %>%
  dplyr::group_by(sample, chain, call,gene_type) %>%
  dplyr::summarise(allele_count = n_distinct(allele))


combined_data <- allele_count_sample %>%
  full_join(allele_count_repertior_true, by = c("sample", "call", "chain", "gene_type"))

combined_data <- combined_data %>%
  dplyr::rename(
    personal_refrence = allele_count, 
    present_in_repertoire = allele_count_true
  )

combined_data_long <- combined_data %>%
  pivot_longer(
    cols = c(personal_refrence, present_in_repertoire),
    names_to = "source",
    values_to = "allele_count_combined",
    names_repair = "unique"  # This repairs any duplicate names
  )


combined_data_long <- combined_data_long %>%
  dplyr::mutate(source = factor(source, levels = c("personal_refrence", "present_in_repertoire")))


plot4 <- ggplot(combined_data_long, aes(x = source, y = allele_count_combined, fill = source)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +  # Use fill aesthetic for coloring
  labs(x = "\nGG\n\n Expressed\n\t(GG)\n", y = "") +  # Remove default x-axis title
  scale_x_discrete(
    labels = c(
      "personal_refrence" = "\tGG",
      "present_in_repertoire" = "\tExpressed\n\t(GG)"
    )
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA)) + 
  facet_wrap(~gene_type, scales = "free_y", ncol = 1) +  # Facet by gene_type
  theme(
    axis.text.y = element_blank(),
    #axis.text.x = element_text(size = 34, angle = 90, vjust = 0.5),          # Adjust x-axis text size for better readability
    axis.text.x = element_blank(), 
    #axis.title = element_text(size = 34, angle = 90),  
    axis.title = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 10, b = 10, l = 0, r = -20),
    plot.title = element_text(size = 24),
    panel.spacing = unit(2, "lines")
  ) +
  scale_fill_manual(
    values = c(
      "personal_refrence" = gg_color, 
      "present_in_repertoire" = personal_color
    )
  )

## plot 4 summary

combined_data_long %>%
  dplyr::group_by(gene_type, sample) %>%
  dplyr::summarise(frac = allele_count_combined[source == "present_in_repertoire"] /
                     allele_count_combined[source != "present_in_repertoire"], .groups = "drop") %>%
  dplyr::group_by(gene_type) %>%
  dplyr::summarise(
    median_frac = median(frac, na.rm = TRUE),
    lower_ci = quantile(frac, 0.025, na.rm = TRUE),
    upper_ci = quantile(frac, 0.975, na.rm = TRUE)
  )


combined_data_long %>%
  dplyr::group_by(gene_type, sample) %>%
  dplyr::summarise(frac =allele_count_combined[source != "present_in_repertoire"], .groups = "drop") %>%
  dplyr::group_by(gene_type) %>%
  dplyr::summarise(
    median = median(frac, na.rm = TRUE),
    lower = min(frac, na.rm = TRUE),
    upper = max(frac, na.rm = TRUE)
  )


combined_data_long %>%
  dplyr::group_by(gene_type, sample) %>%
  dplyr::summarise(frac =allele_count_combined[source == "present_in_repertoire"], .groups = "drop") %>%
  dplyr::group_by(gene_type) %>%
  dplyr::summarise(
    median = median(frac, na.rm = TRUE),
    lower = min(frac, na.rm = TRUE),
    upper = max(frac, na.rm = TRUE)
  )


plot1 <- plot1 + theme(plot.margin = unit(c(0, 4, 0, 0), "pt"))  # Add margin to the right
plot4 <- plot4 + theme(plot.margin = unit(c(0, 0, 0, -50), "pt")) # Add negative margin to the left
plot_cdf <- plot_cdf + theme(plot.margin = unit(c(0, 0, 0, -10), "pt"))


combined_plot <- plot_strip_only + 
  plot1 + 
  plot4 + 
  wrap_elements(plot3)  + 
  free(wrap_elements(plot_cdf)) + 
  plot_layout(nrow = 1, widths = c(0.01,3.3, 0.8,1.8, 2)) +
  plot_annotation(tag_levels = list(c("","A", "", "B", "C"))) +
  inset_element(legend_plot2, left = -1.5,
                bottom = 0.035,
                right = 0,
                top = 0,
                align_to = "plot", ignore_tag = T, on_top = F) &
  theme(plot.tag = element_text(size = 36, hjust = -1))



ggsave(filename = "figures/figure2.pdf", plot = combined_plot, 
       width = 30, height = 25, dpi = 300,bg = "white",limitsize = FALSE)



