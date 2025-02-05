# themes=upset_modify_themes(
#   list(
#     'intersections_matrix'=theme(text=element_text(size=16), 
#                                  axis.title.x = element_blank(),
#                                  axis.text.x = element_blank(),
#                                  axis.ticks.length.x=unit(0,'lines'),
#                                  panel.grid.minor.y = element_blank())
#   )
# )
# 
# base_annotations=list(
#   '# of shared Heptamer-Nonamer'=intersection_size2(
#     counts=TRUE,
#     #mapping=aes(fill=rss_group),
#     text_colors = c("black","black")
#   )+theme(axis.text.y=element_text(size=14),
#           axis.title.y = element_text(size=16),
#           legend.text = element_text(size=20),
#           legend.title = element_text(size=21)) +
#     guides(fill="none") +
#     labs(fill="") + coord_cartesian(clip = "off")
# )


upset_adapted_v2 = function (data, data2, intersect, base_annotations = "auto", name = "group", 
                             annotations = list(), themes = upset_themes, stripes = upset_stripes(), 
                             labeller = identity, height_ratio = 0.5, width_ratio = 0.3, 
                             wrap = FALSE, set_sizes = upset_set_size(), mode = "distinct", 
                             queries = list(), guides = NULL, encode_sets = TRUE, matrix = intersection_matrix(), 
                             ...) 
{
  if (!is.null(guides)) {
    check_argument(guides, allowed = c("keep", "collect", 
                                       "over"), "guides")
  }
  mode = ComplexUpset:::solve_mode(mode)
  if (class(base_annotations) == "character") {
    if (base_annotations != "auto") {
      stop("Unsupported value for `base_annotations`: provide a named list, or `\"auto\"`")
    }
    else {
      base_annotations = list(`Intersection size` = intersection_size(counts = TRUE, 
                                                                      mode = mode))
    }
  }
  if (class(stripes) != "upset_stripes") {
    stripes = upset_stripes(colors = stripes)
  }
  annotations = c(annotations, base_annotations)
  data = upset_data(data, intersect, mode = mode, encode_sets = encode_sets)
  intersections_sorted = rev(data$sorted$intersections)
  intersections_limits = intersections_sorted[intersections_sorted %in% 
                                                data$plot_intersections_subset]
  
  sets_limits = data$sorted$groups[data$sorted$groups %in% 
                                     data$plot_sets_subset]
  show_overall_sizes = !(inherits(set_sizes, "logical") && 
                           set_sizes == FALSE)
  matrix_intersect_queries = ComplexUpset:::intersect_queries(ComplexUpset:::queries_for(queries, 
                                                           "intersections_matrix"), data)
  matrix_group_by_queries = ComplexUpset:::group_by_queries(ComplexUpset:::queries_for(queries, 
                                                         "intersections_matrix"), data$sanitized_labels)
  matrix_set_queries = ComplexUpset:::set_queries(ComplexUpset:::queries_for(queries, "intersections_matrix"), 
                                   data$sanitized_labels)
  intersection_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 
                                                  "intersection", matrix_intersect_queries)
  group_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 
                                           "group_by_group", matrix_group_by_queries)
  set_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 
                                         "group", matrix_set_queries)
  query_matrix = ComplexUpset:::merge_rows(ComplexUpset:::merge_rows(intersection_query_matrix, 
                                       group_query_matrix), set_query_matrix)
  query_matrix = query_matrix[query_matrix$value == TRUE, ]
  matrix_frame = data$matrix_frame[data$matrix_frame$group %in% 
                                     data$plot_sets_subset, ]
  
  ### we need to create N columns based on the max number of shared couples plus 1.
  ## we then need to fill the proportion of the alleles in each category
  ## the rest should be zeros
  ## the size of the pie will be the log number of alleles or bin number of alleles.
  ## the colors of the upper bar will reflect the category.
  ## the plus one will reflect situation where there is no shared alleles. 
  ## maybe we can just do them as points. 
  unique_intersect <- unique(grep("-",matrix_frame$intersection,value=T))
  max_number_of_groups <- max(data$sizes$exclusive_intersection[unique_intersect])
  groups_to_add <- max_number_of_groups + 1
  ## add groups_to_add as columns to the matrix_frame
  matrix_frame_copy <- matrix_frame
  for (i in 1:groups_to_add){
    matrix_frame_copy[[paste0("group_",i)]] <- 0
  }
  
  ## group_1 should be filed just with the values of singles
  tmp <- matrix_frame_copy[!grepl("-",matrix_frame_copy$intersection) & matrix_frame_copy$value==T,]
  for(i in rownames(tmp)){
    g <- tmp[i,"intersection"]
    rss_sequences <- data$with_sizes[data$with_sizes$intersection==g,"rss_sequence"]
    matrix_frame_copy[i,"group_1"] <- data2[rss_sequence %in% rss_sequences & gene_type==data$non_sanitized_labels[as.character(g)],sum(count)]
  }
  
  ## now do the same for the multi intersections. 
  rss_sequences_groups <- list()
  for(g in unique_intersect){
    rss_sequences <- unique(data$with_sizes[data$with_sizes$exclusive_intersection==g,"rss_sequence"])
    ii <- 2
    for(s in rss_sequences){
      tmp_count <- data2[rss_sequence == s & gene_type %in% data$non_sanitized_labels[unlist(strsplit(as.character(g),"-"))],]
      for(g_i in 1:nrow(tmp_count)){
        g_type <- tmp_count$gene_type[g_i]
        group <- data$sanitized_labels[g_type]
        matrix_frame_copy[matrix_frame_copy$value==T & matrix_frame_copy$group==group & matrix_frame_copy$intersection==g,paste0("group_",ii)] <- tmp_count$count[g_i]
      }
      rss_sequences_groups[[s]] <- paste0("group_",ii)
      ii <- ii+1
    }
  }
  
  group_cols <- paste0("group_", 1:groups_to_add)
  
  df_long <- matrix_frame_copy %>%
    pivot_longer(
      cols = all_of(group_cols),
      names_to = "pie",
      values_to = "count"
    )
  
  df_long$count <- as.numeric(df_long$count)
  
  df_long <- df_long %>%
    group_by(intersection, group) %>%
    mutate(
      total = sum(count, na.rm = TRUE),
      fraction = count / total
    ) %>%
    ungroup()
  
  df_long$fraction[is.na(df_long$fraction)] <- 0
  
  df_long <- df_long %>%
    group_by(intersection, group) %>%
    arrange(pie) %>%
    mutate(
      end = if_else(value, 2 * pi * cumsum(fraction), NA_real_),
      start = if_else(value, lag(end, default = 0), NA_real_)
    ) %>%
    ungroup()
  
  
  df_long <- df_long %>%
    mutate(
      intersection_or = intersection,
      group_or = group,
    )
  
  df_long$tags <- data$non_sanitized_labels[df_long$group_or]
  
  df_long <- df_long %>%
    mutate(
      intersection = match(as.character(intersection), intersections_limits),
      group = match(as.character(group), sets_limits)
    )
  
  df_long$total_color <- df_long$total
  df_long$total_color[df_long$total_color==0] <- NA
  # breaks <- c(0, 10, 50, 100, 200, 400, 600, Inf)
  # labels <- c("1-10", "11-50", "51-100", "101-200", "201-400", "401-600", "600+")
  # df_long$total_color_bins <- cut(
  #   df_long$total_color,
  #   breaks = breaks,
  #   labels = labels,
  #   right = TRUE,
  #   include.lowest = TRUE
  # )
  
  breaks <- c(1, 2, 11, 101, 200, 400, 800)
  #labels <- c("1", "2-10", "11-100", "101-200", "201-400", "401-800")
  
  df_long$total_color_bins <- cut(
    df_long$total_color,
    breaks = breaks,
    #labels = labels,
    right = TRUE,
    include.lowest = TRUE
  )
  
  matrix = intersection_matrix_update()
  intersections_matrix = matrix %+% df_long
  point_geom = intersections_matrix$geom
  if (!is.null(point_geom$aes_params$size)) {
    dot_size = point_geom$aes_params$size
  }else {
    dot_size = 1
  }
  
  geom_layers = c(
    #interconnectors on the dots
    list(
      ComplexUpset:::`*.gg`(intersections_matrix$segment, geom_segment(
        aes(
          x = intersection,
          xend = intersection,
          y = segment_end_long(df_long, intersection, min),  # Replace `head`
          yend = segment_end_long(df_long, intersection, max) # Replace `tail`
        ),
        na.rm = TRUE
      ))
    ),
    # the dots outline
    list(
      ComplexUpset:::`*.gg`(intersections_matrix$geom, geom_point(
        mapping = aes(x = intersection, y = group, color = (total_color_bins)),
        # color=ifelse(
        #   df_long$value,
        #   intersections_matrix$outline_color$active,
        #   intersections_matrix$outline_color$inactive
        # ),
        size = 10,
        na.rm = TRUE
      ))
    ),
    
    # Replace the dot with pie charts
    list(
      ComplexUpset:::`*.gg`(
        ggforce::geom_arc_bar(),
        ggforce::geom_arc_bar(
          data = df_long[df_long$value==TRUE,],
          aes(
            x0 = intersection,
            y0 = group,
            r0 = 0.2,
            r = 0.3,
            start = start,
            end = end,
            fill = pie
          ),
          na.rm = TRUE,
          inherit.aes = F,
          color = "transparent", show.legend = F
        )
      )
    )
  )
  
  intersections_matrix$layers = c(matrix_background_stripes_long(df_long, stripes), geom_layers, intersections_matrix$layers)
  
  
  
  scale_intersections <- scale_x_continuous(
    breaks = seq_along(intersections_limits),  # Numeric positions
    labels = intersections_limits    # Custom labels
  )
  
  set_labels <- unique(data$non_sanitized_labels[sets_limits])
  
  y_scale <- scale_y_continuous(
    breaks = seq_along(sets_limits),            # Numeric positions
    labels = gsub("IGHD_5","5'IGHD",gsub("IGHD_3","3'IGHD",set_labels))              # Custom labels
  )
  
  user_y_scale = ComplexUpset:::get_scale(intersections_matrix, "y")
  if (length(user_y_scale) == 0) {
    user_y_scale = scale_y_discrete()
  }else {
    user_y_scale = user_y_scale[[1]]
    user_y_scale$limits = y_scale$limits
    user_y_scale$labels = y_scale$labels
    y_scale = NULL
  }
  matrix_default_colors = list(`TRUE` = "black", `FALSE` = "grey85")
  matrix_guide = "none"
  matrix_breaks = names(matrix_default_colors)
  if (!is.null(names(stripes$colors))) {
    matrix_default_colors = c(matrix_default_colors, stripes$colors)
    matrix_guide = guide_legend()
    matrix_breaks = names(stripes$colors)
  }
  
  colors <- hcl.colors(length(levels(df_long$total_color_bins)), palette = "Heat 2", rev = T)
  
  pl <- setNames(hcl.colors(5, palette = "Zissou1"),unique(df_long$pie))
  
  intersections_matrix = (intersections_matrix + xlab(name) + 
                            scale_intersections + y_scale + 
                            scale_color_manual(
                              values = setNames(colors,levels(df_long$total_color_bins)),
                              breaks = levels(df_long$total_color_bins),
                              na.value = intersections_matrix$outline_color$inactive,
                              name = "# Alleles"
                            ) + scale_fill_manual(
                              values = pl
                            ) )
  
  # get the legends
  leg <- ggpubr::get_legend(intersections_matrix + 
                              theme(
                                legend.text = element_text(size = 30),
                                legend.title = element_text(size = 30, vjust = 1, hjust = 1),
                                legend.direction = "horizontal",
                                legend.background = element_blank(),
                                legend.key = element_rect(fill = "transparent", linewidth = 1, colour = "gray90")
                              ) +
                              guides(fill="none", color = guide_legend(nrow = 3)) #, color =  guide_colourbar(barwidth = 15, barheight = 2)
                            )
  leg <- ggpubr::as_ggplot(leg)
  
  intersections_matrix = (intersections_matrix + themes$intersections_matrix + guides(fill="none", color = "none"))
  rows = list()
  if (show_overall_sizes) {
    is_set_size_on_the_right = !is.null(set_sizes$position) && 
      set_sizes$position == "right"
  }
  annotation_number = 1
  for (name in names(annotations)) {
    annotation = annotations[[name]] + scale_fill_manual(values=pl)
    geoms = annotation$geom
    annotation_mode = mode
    for (layer in annotation$layers) {
      if (inherits(layer$stat, "StatMode")) {
        annotation_mode = layer$stat_params$mode
      }
    }
    annotation_data = data$with_sizes[data$with_sizes[ComplexUpset:::get_mode_presence(annotation_mode, 
                                                                        symbol = FALSE)] == 1, ]
    annotation_data$rss_group <- sapply(annotation_data$rss_sequence, function(s){
      if(s %in% names(rss_sequences_groups)){
        rss_sequences_groups[[s]]
      }else{
        NA_character_
      }
    })
    ## add the contiues scale
    
    annotation_data$intersection = match(as.character(annotation_data$intersection), intersections_limits)
    
    if (!inherits(geoms, "list")) {
      geoms = list(geoms)
    }
    annotation_queries = ComplexUpset:::intersect_queries(ComplexUpset:::queries_for(queries, 
                                                       name), data)
    if (nrow(annotation_queries) != 0) {
      highlight_data = merge(annotation_data, annotation_queries, 
                             by.x = "intersection", by.y = "intersect", all.y = TRUE)
      if (is.null(annotation$highlight_geom)) {
        highlight_geom = geoms
      }
      else {
        highlight_geom = annotation$highlight_geom
        if (!inherits(highlight_geom, "list")) {
          highlight_geom = list(highlight_geom)
        }
      }
      geoms_plus_highlights = ComplexUpset:::add_highlights_to_geoms(geoms, 
                                                      highlight_geom, highlight_data, annotation_queries)
    }else {
      geoms_plus_highlights = geoms
    }
    if (!is.null(annotation$top_geom)) {
      geoms_plus_highlights = c(geoms_plus_highlights, 
                                annotation$top_geom)
    }
    if (name %in% names(themes)) {
      selected_theme = themes[[name]]
    }else {
      selected_theme = themes[["default"]]
    }
    if (!is.null(guides) && guides == "over" && ceiling(length(annotations)/2) == 
        annotation_number) {
      spacer = guide_area()
    }else {
      spacer = plot_spacer()
    }
    if (show_overall_sizes && !is_set_size_on_the_right) {
      rows[[length(rows) + 1]] = spacer
    }
    if (is.ggplot(annotation)) {
      if (is.null(annotation$mapping$x)) {
        annotation = annotation + aes(x = intersection)
      }
      annotation_plot = annotation %+% annotation_data
      user_theme = annotation_plot$theme
      annotation_plot = annotation_plot + selected_theme + 
        do.call(theme, user_theme)
      if (is.null(annotation_plot$default_y) && !is.null(annotation_plot$mapping$y)) {
        annotation_plot$default_y = quo_name(annotation_plot$mapping$y)
      }
      if (is.null(annotation_plot$labels$y) || (!is.null(annotation_plot$default_y) && 
                                                annotation_plot$default_y == annotation_plot$labels$y)) {
        annotation_plot = annotation_plot + ylab(name)
      }
    }else {
      annotation_plot = ggplot(annotation_data, annotation$aes) + 
        selected_theme + xlab(name) + ylab(name)
    }
    user_layers = annotation_plot$layers
    annotation_plot$layers = c()
    annotation_plot = annotation_plot + geoms_plus_highlights
    annotation_plot$layers = c(annotation_plot$layers, user_layers)
    rows[[length(rows) + 1]] = (annotation_plot + scale_intersections)
    if (show_overall_sizes && is_set_size_on_the_right) {
      rows[[length(rows) + 1]] = spacer
    }
    annotation_number = annotation_number + 1
  }
  if (show_overall_sizes) {
    set_sizes_data = data$presence[data$presence$group %in% 
                                     data$plot_sets_subset, ]
    if (set_sizes$filter_intersections) {
      set_sizes_data = set_sizes_data[set_sizes_data$intersection %in% 
                                        data$plot_intersections_subset, ]
    }
    
    set_sizes_data$group = match(as.character(set_sizes_data$group), sets_limits)
    
    overall_sizes_queries = ComplexUpset:::set_queries(ComplexUpset:::queries_for(queries, 
                                                    "overall_sizes"), data$sanitized_labels)
    overall_sizes_highlights_data = ComplexUpset:::get_highlights_data(set_sizes_data, 
                                                        "group", overall_sizes_queries)
    if (nrow(overall_sizes_queries) != 0) {
      highlight_geom = set_sizes$highlight_geom
      if (!inherits(highlight_geom, "list")) {
        highlight_geom = list(highlight_geom)
      }
      overall_sizes_highlights_data$group = factor(overall_sizes_highlights_data$group)
      geom = ComplexUpset:::add_highlights_to_geoms(set_sizes$geom, highlight_geom, 
                                     overall_sizes_highlights_data, overall_sizes_queries, 
                                     kind = "group")
    }else {
      geom = set_sizes$geom
    }
    if (is_set_size_on_the_right) {
      default_scale = scale_y_continuous()
    }else {
      default_scale = scale_y_reverse()
    }
    set_sizes$layers = c(ComplexUpset:::matrix_background_stripes(data, 
                                                   stripes, "vertical"), geom, set_sizes$layers)
    overall_sizes = (set_sizes %+% set_sizes_data + aes(x = group) + 
                       themes$overall_sizes + do.call(theme, set_sizes$theme) + 
                       coord_flip() + scale_x_discrete(limits = sets_limits) + 
                       ComplexUpset:::scale_if_missing(set_sizes, axis = "y", scale = default_scale) + 
                       ComplexUpset:::scale_if_missing(set_sizes, "colour", scale_color_manual(values = matrix_default_colors, 
                                                                                guide = "none")))
    if (is_set_size_on_the_right) {
      matrix_row = list(intersections_matrix, overall_sizes)
    }else {
      matrix_row = list(overall_sizes, intersections_matrix)
    }
  }else {
    matrix_row = list(intersections_matrix)
  }
  if (length(rows)) {
    annotations_plots = Reduce(f = "+", rows)
    matrix_row = c(list(annotations_plots), matrix_row)
  }else {
    annotations_plots = list()
  }
  plot = Reduce(f = "+", matrix_row)
  if (show_overall_sizes) {
    if (is_set_size_on_the_right) {
      width_ratio = 1 - width_ratio
    }
    width_ratios = c(width_ratio, 1 - width_ratio)
  }else {
    width_ratios = 1
  }
  if (!is.null(guides) && guides == "over") {
    guides = "collect"
  }
  plot = plot + plot_layout(widths = width_ratios, ncol = 1 + 
                              ifelse(show_overall_sizes, 1, 0), nrow = length(annotations) + 
                              1, heights = c(rep(1, length(annotations)), height_ratio), 
                            guides = guides)
  
  if (wrap) {

    return(list(plot = wrap_elements(plot), leg = leg))
  }
  else {
    return(list(plot = plot, leg = leg))
  }
}


segment_end_long <- function(data_long, intersection, end) {
  # Filter rows for the given intersection
  corresponding <- data_long[data_long$intersection == intersection, ]
  
  # Ensure `group` and `end` are numeric
  corresponding <- corresponding %>%
    mutate(
      group = as.numeric(group),
      fraction = as.numeric(fraction)
    )
  
  # Determine the segment ends
  h <- sapply(seq_len(nrow(corresponding)), function(i) {
    row <- corresponding[i, ]
    group_members <- corresponding %>%
      filter(intersection == row$intersection, count > 0) %>%
      pull(group)
    
    # Use the `end` function (e.g., min, max) on valid numeric group_members
    if (length(group_members) > 0) {
      end(group_members)
    } else {
      NA  # No valid groups
    }
  })
  
  return(h)
}

matrix_background_stripes_long <- function(df_long, stripes, orient = "horizontal") {
  if (!(orient %in% c("horizontal", "vertical"))) {
    stop("Incorrect orient")
  }
  
  # Define aesthetics based on orientation
  if (orient == "horizontal") {
    aes <- aes(x = -Inf, xend = Inf, y = group, yend = group)
  } else {
    aes <- aes(y = -Inf, yend = Inf, x = group, xend = group)
  }
  
  # Get unique groups from df_long
  groups <- sort(unique(df_long$group))
  data <- data.frame(group = groups, group_name = as.character(groups))
  
  # Merge with additional stripe data if provided
  if (!is.null(stripes$data)) {
    data <- merge(data, stripes$data, by.x = "group_name", by.y = "set", all.x = TRUE)
  }
  
  # Create parameters for geom_segment
  params <- list(data = data, mapping = modifyList(aes, stripes$mapping), inherit.aes = F)
  
  # Assign colors if provided
  if (!is.null(stripes$colors) && is.null(names(stripes$colors))) {
    params$color <- rep_len(stripes$colors, length(groups))
  }
  
  # Generate the geom_segment layers
  list(ComplexUpset:::`*.gg`(stripes$geom, do.call(geom_segment, params)))
}

intersection_matrix_update <- function (geom = geom_point(size = 3), segment = geom_segment(), 
                                        outline_color = list(active = "black", inactive = "grey70")) 
{
  plot = ggplot()
  plot$geom = geom
  plot$segment = segment
  plot$outline_color = outline_color
  plot
}


intersection_size2 <- function (mapping = aes(), counts = TRUE, bar_number_threshold = 0.85, 
                                text_colors = c(on_background = "black", on_bar = "white"), 
                                text = list(), text_mapping = aes(), mode = "distinct", 
                                position = position_stack(), ...) 
{
  size = get_size_mode(mode)
  lab = switch(mode, exclusive_intersection = "Intersection size", 
               inclusive_intersection = "Inclusive intersection size", 
               inclusive_union = "Union size", exclusive_union = "Exclusive union size")
  if (counts) {
    text = modifyList(ComplexUpset:::intersection_size_text, text)
    text_mapping = modifyList(aes(label = !!size, y = !!size, colour = ifelse(!!size <= bar_number_threshold * max(!!size, na.rm = TRUE), 
                                                                              "on_background", "on_bar")), text_mapping)
    counts_geoms = list(do.call(geom_text, c(list(stat = "unique", 
                                                  text_mapping, na.rm = TRUE), text)), scale_color_manual(values = text_colors, 
                                                                                                          guide = "none"))
  }
  else {
    counts_geoms = list()
  }
  bar_geom = list(stat_summary(fun = sum, geom = "bar", position = position, 
                               na.rm = TRUE, show.legend=T, ...))
  ComplexUpset:::convert_annotation(aes = modifyList(aes(x = intersection, 
                                                         y = !!ComplexUpset:::get_mode_presence(mode)), mapping), geom = bar_geom, 
                                    highlight_geom = bar_geom, top_geom = counts_geoms) + 
    ylab(lab) + upset_mode(mode)
}



