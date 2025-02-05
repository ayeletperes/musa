
filtered_results <- readRDS("data/rhesus_macaque_data/04_02_d_allele_filtered_results.rds")


allele_medians <- filtered_results %>%
  group_by(allele) %>%
  summarize(
    median_value = median((d_germline_end - d_germline_start + 1) / nchar, na.rm = TRUE)
  ) %>%
  arrange(median_value)

# Create labels for alleles with median_value < 0.5
low_median_alleles <- allele_medians %>%
  filter(median_value <= 0.5 & grepl("_", allele)) %>%
  mutate(allele_label = paste0("D", row_number()))

# Map the sequential labels back to the full dataset
allele_medians <- allele_medians %>%
  left_join(
    low_median_alleles %>% select(allele, allele_label),
    by = "allele"
  ) %>%
  mutate(
    allele_label = if_else(is.na(allele_label), allele, allele_label)
  )

# Add the new labels to filtered_results
filtered_results <- filtered_results %>%
  left_join(allele_medians %>% select(allele, allele_label), by = "allele")


filtered_results <- filtered_results %>%
  mutate(
    group = ifelse(grepl("^D", allele_label), "D* Alleles", "Other Alleles")
  )

p_more_0 <- ggplot(filtered_results, aes(x = allele_label, y = ((d_germline_end - d_germline_start + 1) / nchar), fill = has_underscore)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") + 
  labs(
    x = "Allele",
    y = "Sequence Length / Germline Length",
    fill = " "
  ) +
  scale_fill_manual(
    values = c("novel" = "orange", "baseline" = "blue"),
    guide = guide_legend(
      override.aes = list(fill = c("blue", "orange"), color = NA, linetype = 0)
    )
  ) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 30),
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.key = element_blank()
  ) +
  theme(
    plot.margin = margin(t = 10, r = 20, b = 10, l = 40),
    axis.text.y = element_text(margin = margin(r = 10))
  ) +
  facet_wrap(~group, scales = "free_x", ncol = 1)

d_star_table <- allele_medians %>%
  filter(grepl("^D", allele_label)) %>%
  select(allele_label, allele) %>%
  mutate(index = row_number())  # Add a row index

colnames(d_star_table)

colnames(d_star_table) <- c("Label name", "Allele", "index")

n_rows <- ceiling(nrow(d_star_table) / 3)

# Pad the table with NAs to ensure equal-sized groups
padded_table <- d_star_table %>%
  bind_rows(tibble(
    `Label name` = rep(NA_character_, 3 * n_rows - nrow(d_star_table)),
    Allele = rep(NA_character_, 3 * n_rows - nrow(d_star_table)),
    index = rep(NA_integer_, 3 * n_rows - nrow(d_star_table))
  ))

# Split the data into three equal parts
split_tables <- split(padded_table, (seq_len(nrow(padded_table)) - 1) %/% n_rows)

# Combine the splits side by side
result_table <- bind_cols(
  split_tables[[1]],
  split_tables[[2]],
  split_tables[[3]]
)

# Rename columns to remove numbers
colnames(result_table) <- c("Label name", "Allele","index_1",
                            "Label name", "Allele","index_2",
                            "Label name", "Allele", "index_4")
result_table <- result_table[,c(1,2,4,5,7,8)]

result_table[is.na(result_table)] <- ""

library(gridExtra)
# Prepare the final table for display
d_star_grob <- tableGrob(result_table, rows = NULL)

library(cowplot)

# Combine the plot and table
final_plot <- plot_grid(
  p_more_0,                        # Main ggplot
  ggdraw() + draw_grob(d_star_grob),                     # Table grob
  ncol = 1,                        # Stack vertically
  rel_heights = c(1, 1)            # Adjust relative heights
)


ggsave("figures/04_02_d_filter_plot.pdf", final_plot, width = 30, height = 30, dpi = 300,bg = "white",limitsize = FALSE)

