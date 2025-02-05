upset_adapted = function(
    data,
    base_annotations='auto',
    name='group',
    annotations=list(),
    themes=upset_themes,
    stripes=upset_stripes(),
    labeller=identity,
    height_ratio=0.5,
    width_ratio=0.3,
    wrap=FALSE,
    set_sizes=upset_set_size(),
    mode='distinct',
    queries=list(),
    guides=NULL,
    encode_sets=TRUE,
    matrix=intersection_matrix(),
    add=TRUE,
    limits = c(-600,600),
    ...
) {
  if (!is.null(guides)) {
    check_argument(guides, allowed = c('keep', 'collect', 'over'), 'guides')
  }
  
  mode = ComplexUpset:::solve_mode(mode)
  
  if (inherits(base_annotations, 'character')) {
    if (base_annotations != 'auto') {
      stop('Unsupported value for `base_annotations`: provide a named list, or `"auto"`')
    } else {
      base_annotations = list(
        'Intersection size'=intersection_size(counts=TRUE, mode=mode)
      )
    }
  }
  
  # for backwards compatibility pre 1.2
  if (!inherits(stripes, 'upset_stripes')) {
    stripes = upset_stripes(colors=stripes)
  }
  
  annotations = c(annotations, base_annotations)
  
  intersections_sorted = rev(data$sorted$intersections)
  intersections_limits = intersections_sorted[intersections_sorted %in% data$plot_intersections_subset]
  
  scale_intersections = scale_x_discrete(limits=intersections_limits)
  
  sets_limits = data$sorted$groups[data$sorted$groups %in% data$plot_sets_subset]
  
  show_overall_sizes = !(inherits(set_sizes, 'logical') && set_sizes == FALSE)
  
  matrix_intersect_queries = ComplexUpset:::intersect_queries(ComplexUpset:::queries_for(queries, 'intersections_matrix'), data)
  matrix_group_by_queries = ComplexUpset:::group_by_queries(ComplexUpset:::queries_for(queries, 'intersections_matrix'), data$sanitized_labels)
  matrix_set_queries = ComplexUpset:::set_queries(ComplexUpset:::queries_for(queries, 'intersections_matrix'), data$sanitized_labels)
  
  intersection_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 'intersection', matrix_intersect_queries)
  group_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 'group_by_group', matrix_group_by_queries)
  set_query_matrix = ComplexUpset:::get_highlights_data(data$matrix_frame, 'group', matrix_set_queries)
  
  
  query_matrix = ComplexUpset:::merge_rows(
    ComplexUpset:::merge_rows(
      intersection_query_matrix,
      group_query_matrix
    ),
    set_query_matrix
  )
  
  query_matrix = query_matrix[query_matrix$value == TRUE, ]
  
  matrix_frame = data$matrix_frame[data$matrix_frame$group %in% data$plot_sets_subset, ]
  intersections_matrix = matrix %+% matrix_frame
  
  point_geom = intersections_matrix$geom
  
  if (!is.null(point_geom$aes_params$size)) {
    dot_size = point_geom$aes_params$size
  } else {
    dot_size = 1
  }
  
  geom_layers = c(
    # the dots outline
    list(
      ComplexUpset:::`*.gg`(intersections_matrix$geom, geom_point(
        color=ifelse(
          matrix_frame$value,
          intersections_matrix$outline_color$active,
          intersections_matrix$outline_color$inactive
        ),
        size=dot_size * 7/6,
        na.rm=TRUE
      ))
    ),
    # the dot
    list(
      ComplexUpset:::`*.gg`(intersections_matrix$geom, geom_point(
        aes(color=value),
        size=dot_size,
        na.rm=TRUE
      ))),
    # interconnectors on the dots
    list(
      ComplexUpset:::`*.gg`(intersections_matrix$segment, geom_segment(aes(
        x=intersection,
        xend=intersection,
        y=ComplexUpset:::segment_end(matrix_frame, data, intersection, head),
        yend=ComplexUpset:::segment_end(matrix_frame, data, intersection, tail)
      ), na.rm=TRUE))
    ),
    # highlighted interconnectors
    ComplexUpset:::highlight_layer(
      intersections_matrix$segment,
      geom_segment,
      query_matrix,
      args=list(
        aes(
          x=intersection,
          xend=intersection,
          y=ComplexUpset:::segment_end(query_matrix, data, intersection, head),
          yend=ComplexUpset:::segment_end(query_matrix, data, intersection, tail)
        ),
        na.rm=TRUE
      )
    ),
    # the highlighted dot
    ComplexUpset:::highlight_layer(
      intersections_matrix$geom,
      geom_point,
      query_matrix,
      args=list(size=dot_size)
    )
  )
  
  intersections_matrix$layers = c(
    ComplexUpset:::matrix_background_stripes(data, stripes),
    geom_layers,
    intersections_matrix$layers
  )
  y_scale = scale_y_discrete(
    limits=sets_limits,
    labels=function(sets) { labeller(data$non_sanitized_labels[sets]) }
  )
  
  user_y_scale = ComplexUpset:::get_scale(intersections_matrix, 'y')
  
  if (length(user_y_scale) == 0) {
    user_y_scale = scale_y_discrete()
  } else {
    user_y_scale = user_y_scale[[1]]
    user_y_scale$limits = y_scale$limits
    user_y_scale$labels = y_scale$labels
    y_scale = NULL
  }
  
  matrix_default_colors = list('TRUE'='black', 'FALSE'='grey85')
  matrix_guide = "none"
  matrix_breaks = names(matrix_default_colors)
  if (!is.null(names(stripes$colors))) {
    matrix_default_colors = c(
      matrix_default_colors,
      stripes$colors
    )
    matrix_guide = guide_legend()
    matrix_breaks = names(stripes$colors)
  }
  
  intersections_matrix = (
    intersections_matrix
    + xlab(name)
    + scale_intersections
    + y_scale
    + ComplexUpset:::scale_if_missing(
      intersections_matrix,
      'colour',
      scale_color_manual(
        values=matrix_default_colors,
        guide=matrix_guide,
        breaks=matrix_breaks
      )
    )
    + themes$intersections_matrix
  )
  
  rows = list()
  
  if (show_overall_sizes) {
    is_set_size_on_the_right = !is.null(set_sizes$position) && set_sizes$position == 'right'
  }
  
  annotation_number = 1
  
  for (name in names(annotations)) {
    annotation = annotations[[name]]
    
    geoms = annotation$geom
    
    annotation_mode = mode
    
    for (layer in annotation$layers) {
      if (inherits(layer$stat, 'StatMode')) {
        annotation_mode = layer$stat_params$mode
      }
    }
    
    annotation_data = data$with_sizes[data$with_sizes[ComplexUpset:::get_mode_presence(annotation_mode, symbol=FALSE)] == 1, ]
    
    if (!inherits(geoms, 'list')) {
      geoms = list(geoms)
    }
    
    annotation_queries = ComplexUpset:::intersect_queries(ComplexUpset:::queries_for(queries, name), data)
    
    if (nrow(annotation_queries) != 0) {
      highlight_data = merge(annotation_data, annotation_queries, by.x='intersection', by.y='intersect', all.y=TRUE)
      
      if (is.null(annotation$highlight_geom)) {
        highlight_geom = geoms
      } else {
        highlight_geom = annotation$highlight_geom
        if (!inherits(highlight_geom, 'list')) {
          highlight_geom = list(highlight_geom)
        }
      }
      
      geoms_plus_highlights = add_highlights_to_geoms(geoms, highlight_geom, highlight_data, annotation_queries)
    } else {
      geoms_plus_highlights = geoms
    }
    
    if (!is.null(annotation$top_geom)) {
      geoms_plus_highlights = c(geoms_plus_highlights, annotation$top_geom)
    }
    
    if (name %in% names(themes)) {
      selected_theme = themes[[name]]
    } else {
      selected_theme = themes[['default']]
    }
    
    if (!is.null(guides) && guides == 'over' && ceiling(length(annotations) / 2) == annotation_number) {
      spacer = guide_area()
    } else {
      spacer = plot_spacer()
    }
    
    if (show_overall_sizes && !is_set_size_on_the_right) {
      rows[[length(rows) + 1]] = spacer
    }
    
    if (is.ggplot(annotation)) {
      if (is.null(annotation$mapping$x)) {
        annotation = annotation + aes(x=intersection)
      }
      annotation_plot = annotation %+% annotation_data
      user_theme = annotation_plot$theme
      annotation_plot = annotation_plot + selected_theme + do.call(theme, user_theme)
      
      if (is.null(annotation_plot$default_y) && !is.null(annotation_plot$mapping$y)) {
        annotation_plot$default_y = quo_name(annotation_plot$mapping$y)
      }
      if (
        is.null(annotation_plot$labels$y)
        ||
        (
          !is.null(annotation_plot$default_y)
          &&
          annotation_plot$default_y == annotation_plot$labels$y
        )
      ) {
        annotation_plot = annotation_plot + ylab(name)
      }
    } else {
      annotation_plot = ggplot(annotation_data, annotation$aes) + selected_theme + xlab(name) + ylab(name)
    }
    
    user_layers = annotation_plot$layers
    annotation_plot$layers = c()
    annotation_plot = annotation_plot + geoms_plus_highlights
    annotation_plot$layers = c(annotation_plot$layers, user_layers)
    
    rows[[length(rows) + 1]] = (
      annotation_plot
      + scale_intersections
    )
    
    if (show_overall_sizes && is_set_size_on_the_right) {
      rows[[length(rows) + 1]] = spacer
    }
    
    annotation_number =  annotation_number + 1
  }
  
  if (show_overall_sizes) {
    set_sizes_data = data$presence[data$presence$group %in% data$plot_sets_subset, ]
    
    if (set_sizes$filter_intersections) {
      set_sizes_data = set_sizes_data[set_sizes_data$intersection %in% data$plot_intersections_subset, ]
    }
    
    overall_sizes_queries = ComplexUpset:::set_queries(ComplexUpset:::queries_for(queries, 'overall_sizes'), data$sanitized_labels)
    overall_sizes_highlights_data = ComplexUpset:::get_highlights_data(set_sizes_data, 'group', overall_sizes_queries)
    if (nrow(overall_sizes_queries) != 0) {
      highlight_geom = set_sizes$highlight_geom
      if (!inherits(highlight_geom, 'list')) {
        highlight_geom = list(highlight_geom)
      }
      overall_sizes_highlights_data$group = factor(overall_sizes_highlights_data$group)
      geom = add_highlights_to_geoms(
        set_sizes$geom,
        highlight_geom,
        overall_sizes_highlights_data,
        overall_sizes_queries,
        kind='group'
      )
    } else {
      geom = set_sizes$geom
    }
    
    if (is_set_size_on_the_right) {
      default_scale = scale_y_continuous()
    } else {
      default_scale = scale_y_reverse()
    }
    
    convert_annotation = function(...) {
      arguments = list(...)
      if (is.null(arguments$aes)) {
        intersection_plot = ggplot()
      } else {
        intersection_plot = ggplot(mapping=arguments$aes)
      }
      intersection_plot$highlight_geom = arguments$highlight_geom
      intersection_plot$top_geom = arguments$top_geom
      intersection_plot$geom = arguments$geom
      intersection_plot
    }
    
    
    set_sizes$layers = c(
      #ComplexUpset:::matrix_background_stripes(data, stripes, 'vertical'),
      geom,
      set_sizes$layers
    )
    
    set_sizes_data_set = melt(table(set_sizes_data$group, set_sizes_data$present))
    names(set_sizes_data_set) = c('group', 'present', 'value')
    set_sizes_data_set$value[set_sizes_data_set$present != 'Not in cohort'] <- set_sizes_data_set$value[set_sizes_data_set$present != 'Not in cohort'] * -1
    set_sizes_data_set$opacity = ifelse(set_sizes_data_set$present!='Not in cohort', 1, 0.5)
    
    x_val <- length(unique(set_sizes_data_set$group))
    yend_left <- max(set_sizes_data_set$value)
    yend_right <- min(set_sizes_data_set$value)
    yend_right_ops <- abs(yend_right)
    
    ## who ever has the smallest absolute will be the factor.
    yend <- min(abs(c(yend_left, yend_right)))
    
    yend_left <- yend
    yend_right <- -yend
    
    if(add){
      overall_sizes = (
        set_sizes %+% set_sizes_data_set
        + aes(x=group, y=value)
        + themes$overall_sizes
        + do.call(theme, set_sizes$theme)
        + coord_flip(clip = "off", xlim = c(1, length(sets_limits)))
        #+ scale_x_discrete(limits=sets_limits)
        + scale_y_continuous(labels = abs, limits = limits)
        + scale_alpha_continuous(range = c(1, 1))
        + guides(fill="none", alpha="none")
        + labs(y="# of Alleles", x = "Source")
        + ComplexUpset:::scale_if_missing(
          set_sizes,
          'colour',
          scale_color_manual(
            values=matrix_default_colors,
            guide="none"
          ))
        + geom_segment(aes(x = x_val+0.5, 
                         xend = x_val+0.5, y = 0, 
                         yend = yend_left * 0.6), 
                     arrow = arrow(length = unit(0.2, "cm"), type="closed"), inherit.aes = FALSE, color="gray80", linewidth=1.5) 
        + geom_segment(aes(x = x_val+0.5, 
                           xend = x_val+0.5, y = 0, 
                           yend = yend_right * 0.6), 
                       arrow = arrow(length = unit(0.2, "cm"), type="closed"), color="gray40", inherit.aes = FALSE, linewidth=1.5)
        + annotate("text", x = x_val+2, 
                   y = yend_left * 0.8, label = "Not in\ncohort", size=5)
        + annotate("text", x = x_val+2, 
                   y = yend_right * 0.8, label = "In\ncohort", size=5)
        + geom_segment(aes(x = 0.7, 
                           xend = x_val+0.3, y = 0, 
                           yend = 0), 
                       inherit.aes = FALSE, color="black", linewidth=0.5) 
      )
    }else{
      overall_sizes = (
        set_sizes %+% set_sizes_data_set
        + aes(x=group, y=value)
        + themes$overall_sizes
        + do.call(theme, set_sizes$theme)
        + coord_flip(clip = "off", xlim = c(1, length(sets_limits)))
        #+ scale_x_discrete(limits=sets_limits)
        + scale_y_continuous(labels = abs, limits = limits)
        + scale_alpha_continuous(range = c(1, 1))
        + guides(fill="none", alpha="none")
        + labs(y="# of Alleles", x = "Source")
        + ComplexUpset:::scale_if_missing(
          set_sizes,
          'colour',
          scale_color_manual(
            values=matrix_default_colors,
            guide="none"
          ))
        + geom_segment(aes(x = 0.7, 
                           xend = x_val+0.3, y = 0, 
                           yend = 0), 
                       inherit.aes = FALSE, color="black", linewidth=0.5) 
      )
      
    }
    
    
    
    if (is_set_size_on_the_right) {
      matrix_row = list(intersections_matrix, overall_sizes)
    } else {
      # on the left by default
      matrix_row = list(overall_sizes, intersections_matrix)
    }
  } else {
    matrix_row = list(intersections_matrix)
  }
  
  ### try to add a venn diaggram
  # arranged = ComplexUpset::arrange_venn(venn_data, sets=c("Inferred", "Broad", "Novel"))
  # 
  # overall_venn = (
  #   ggplot(arranged)
  #   + theme_void()
  #   + coord_fixed(clip = 'off')
  #   + geom_venn_region(data=venn_data, sets = c("Inferred", "Broad", "Novel"), alpha=0.05)
  #   + geom_point(aes(x=x, y=y, color=region, shape=present), size=2.7)
  #   + geom_venn_circle(venn_data, sets=c("Inferred", "Broad", "Novel"), size=1)
  #   + geom_venn_label_set(venn_data, sets=c("Inferred", "Broad", "Novel"), 
  #                         aes(label=region), 
  #                         outwards_adjust=2.6
  #   )
  #   + geom_venn_label_region(venn_data, sets=c("Inferred", "Broad", "Novel"), 
  #                            aes(label=size), 
  #                            outwards_adjust=1.75,
  #                            position=position_nudge(y=0.2))
  #   + scale_color_venn_mix(venn_data, sets=c("Inferred", "Broad", "Novel"), guide='none')
  #   + scale_fill_venn_mix(venn_data[,c("Inferred", "Broad", "Novel")], guide='none')
  #   + scale_shape_manual(
  #     values=c(
  #       'TRUE'='triangle filled',
  #       'FALSE'='triangle down filled'
  #     ),
  #     labels=c(
  #       'TRUE'='Present',
  #       'FALSE'='Not present'
  #     ),
  #     name=''
  #   )
  # )
  # 
  
  if (length(rows)) {
    annotations_plots = Reduce(f='+', rows)
    matrix_row = c(list(annotations_plots), matrix_row)
  } else {
    annotations_plots = list()
  }
  
  plot = Reduce(f='+', matrix_row)
  
  if (show_overall_sizes) {
    if (is_set_size_on_the_right) {
      width_ratio = 1 - width_ratio
    }
    
    width_ratios = c(width_ratio, 1 - width_ratio)
  } else {
    width_ratios = 1
  }
  
  if (!is.null(guides) && guides == 'over') {
    guides = 'collect'  # guide_area() works with collect only
  }
  
  plot = plot + plot_layout(
    widths=width_ratios,
    ncol=1 + ifelse(show_overall_sizes, 1, 0),
    nrow=length(annotations) + 1,
    heights=c(
      rep(1, length(annotations)),
      height_ratio
    ),
    guides=guides
  ) & theme(plot.margin = margin(0, 0, 0, 0, "lines"))
  
  if (wrap) {
    wrap_elements(plot)
  } else {
    plot
  }
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

