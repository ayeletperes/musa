if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  patchwork,
  data.table,
  ComplexUpset,
  ggplot2,
  reshape2
)
# change colors of writing to black. Novel in blue
# seen in repertoires, put the number above. Baseline dark green
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

df_macaque <- fread("data/MUSA_data/MUSA_set_all_2025-02-05_with_asc.csv")

## taking alleles that appeared in the repertoires.

df_macaque_subset <- df_macaque[,.(allele, gene_type, kimdb, `rhgldb+`, imgt, trios, in_baseline_reference,
                                   vrc, guo, novel, sample_count_AIRRseq, sample_count_genomic)]

df_macaque_subset <- df_macaque_subset[,.(
  present_AIRRseq = sum(sample_count_AIRRseq),
  present_Genomic = sum(sample_count_genomic),
  KIMDB = unique(kimdb),
  `RhGLDb+` = unique(`rhgldb+`),
  IMGT = unique(imgt),
  Trios = unique(trios),
  Baseline = any(in_baseline_reference),
  VRC = unique(vrc),
  `Guo et al.` = unique(guo),
  Novel = unique(novel)
), by=.(allele,gene_type)]


df_macaque_subset[,present_AIRRseq:=present_AIRRseq>0]
df_macaque_subset[,present_Genomic:=present_Genomic>0]
df_macaque_subset[,present:=apply(df_macaque_subset[,c("present_AIRRseq","present_Genomic")],1,function(x){
  
  if(x[1]==TRUE & x[2]==TRUE){
    return("Both")
  }else{
    if(x[1]==TRUE){
      return("Repertoire Only (Novel)")
    }else{
      if(x[2]==TRUE){
        return("Genomic Only")
      }else{
        return("Not in cohort")
      }
    }
  }
  
})]
df_macaque_subset$present <- factor(df_macaque_subset$present, levels=c("Both","Genomic Only","Not in cohort","Repertoire Only (Novel)"))


df_macaque_subset[,`:=`(
  present_AIRRseq = NULL,
  present_Genomic = NULL
)]
#df_macaque_subset[,Novel:=grepl("_",allele)]
df_macaque_subset[, (setdiff(names(df_macaque_subset), c("allele", "gene_type","present","Novel","Baseline"))) := 
                    lapply(.SD, function(x) !is.na(x) & x != ""), 
                  .SDcols = setdiff(names(df_macaque_subset), c("allele", "gene_type","present","Novel","Baseline"))
]

df_macaque_subset <- df_macaque_subset[, lapply(.SD, function(x) sum(x)>0), by=.(allele,gene_type,present)]

df_macaque_subset[,gene_type:=substr(allele,1,4)]

## remove any that are only in baseline
other_columns <- setdiff(names(df_macaque_subset), c("allele", "gene_type", "Baseline","present"))

df_macaque_subset[, only_baseline := Baseline & !Reduce(`|`, .SD), .SDcols = other_columns]
df_macaque_subset <- df_macaque_subset[only_baseline==FALSE]
df_macaque_subset <- unique(df_macaque_subset)

## add a column for in database
df_macaque_subset[, in_database := present!="Not in cohort"]

intersect_ = rev(c(
  "Baseline",
  "KIMDB",
  "RhGLDb+",
  "IMGT",
  "Trios",
  "VRC",
  "Guo et al.",
  "Novel"
))


## add a line on zero and then change the opacity of the bar color

## Repertoire only (Novel)

## change the order. 

df_bar <- df_macaque_subset[,.(Baseline=sum(Baseline), Novel=sum(Novel)), by=.(gene_type,present)]
df_bar <- melt.data.table(df_bar, id.vars = c("gene_type","present"), variable.name = "source")
levels <- c("IGHV","IGLV","IGKV","IGHJ","IGLJ","IGKJ","IGHD")
df_bar$gene_type <- factor(df_bar$gene_type, levels=levels)
df_bar$present <- factor(df_bar$present, levels=c("Both","Genomic Only","Not in cohort","Repertoire Only (Novel)"))
df_bar$remove <- FALSE

## edit the IGHD alleles and add a texture for removed

## get total

d_allele_filter <- fread("data/rhesus_macaque_data/04_02_d_allele_more_0_median_value.csv")
in_rep <- d_allele_filter[grepl("_",allele),allele]
in_gg <- df_macaque[gene_type=="IGHD" & novel==TRUE & sample_count_AIRRseq==0 & !(allele %in% in_rep),length(unique(allele))]
filtered_rep <- d_allele_filter[grepl("_",allele) & median_value<=0.5,length(unique(allele))]
# GG - 124
# Both - 127
df_bar <- rbind(df_bar[!(gene_type=="IGHD" & present=="Genomic Only" & source=="Novel")], data.table(
  gene_type = c("IGHD","IGHD"),
  present = c("Genomic Only","Both"),
  source = c("Novel","Novel"),
  value = c(in_gg,filtered_rep),
  remove = c(TRUE, TRUE)
))


library(unikn)
val_bar_fill <- setNames(hcl.colors(4, palette = "Sunset"), c("Both","Repertoire Only (Novel)","Genomic Only","Not in cohort"))

val_bar_fill <- setNames(c("royalblue4","grey70","indianred","royalblue"), c("Both","Genomic Only","Not in cohort","Repertoire Only (Novel)"))

seecol(val_bar_fill)
grays <- rep("gray90", 7)
sets_color <- setNames(c("steelblue1",rev(grays[1:7])),intersect_)
seecol(sets_color)


source("scripts/upset_ComplexUpset.R")

p_upset <- list()
for(gene_type_ in c("IGHV","IGLV","IGKV")){
  df <- df_macaque_subset[gene_type==gene_type_]
  
  ## remove any columns that are all values == FALSE
  df <- df[, lapply(.SD, function(col) if (all(col == FALSE)) NULL else col)]
  
  sets <- intersect(intersect_, colnames(df))
  
  upset_data <- ComplexUpset::upset_data(df, intersect = sets, sort_intersections = F)
  
  ## get number of itersections
  for(source_ in c("IMGT","KIMDB","Trios","RhGLDb+")){
    n1 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection==sprintf("%s-Baseline", source_)& upset_data$with_sizes$in_database==FALSE & upset_data$with_sizes$exclusive_intersection==sprintf("%s-Baseline", source_),])
    n2 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection==sprintf("%s-Baseline", source_)& upset_data$with_sizes$in_database==TRUE & upset_data$with_sizes$exclusive_intersection==sprintf("%s-Baseline", source_),])
    print(sprintf("%s %s not in database: %s", gene_type_, source_, n1))
    print(sprintf("%s %s in database: %s", gene_type_, source_, n2))
    print(sprintf("%s %s total: %s", gene_type_, source_, n2+n1))
  }
  if(gene_type_ == "IGHV"){
    n1 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection=="IMGT-KIMDB-Baseline"& upset_data$with_sizes$in_database==FALSE & upset_data$with_sizes$exclusive_intersection=="IMGT-KIMDB-Baseline",])
    n2 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection=="IMGT-KIMDB-Baseline"& upset_data$with_sizes$in_database==TRUE & upset_data$with_sizes$exclusive_intersection=="IMGT-KIMDB-Baseline",])
    print(sprintf("%s %s not in database: %s", gene_type_, "IMGT-KIMDB", n1))
    print(sprintf("%s %s in database: %s", gene_type_, "IMGT-KIMDB", n2))
    print(sprintf("%s %s total: %s", gene_type_, "IMGT-KIMDB", n2+n1))
  }
  print("----------------------------------")

  
  
  order_ <- sort(upset_data$sizes$exclusive_intersection[!grepl("Novel$",names(upset_data$sizes$exclusive_intersection))], decreasing = F)
  upset_data$sorted$intersections <- c("Novel",names(order_))
  upset_data$sorted$groups <- sets
  
  
  g <- guide_legend(title="")
  g <- if(gene_type_!="IGHV") "none"
  
  overall_sizes_x_label <- if(gene_type_=="IGKV") element_text(size=16)  else element_blank()
  overall_sizes_x_text <- if(gene_type_=="IGKV") element_text(size=16)  else element_blank()
  add <- if(gene_type_=="IGHV") TRUE else FALSE
  p_upset[[gene_type_]] <- upset_adapted(
    data = upset_data,
    base_annotations=list(
      '# of Alleles'=intersection_size2(
        counts=TRUE,
        mapping=aes(fill=present),
        text_colors = c("black","black")
      )+theme(axis.text.y=element_text(size=14),
              axis.title.y = element_text(size=16),
              legend.text = element_text(size=20),
              legend.title = element_text(size=21),
              legend.position = c(0.1, 1.3), 
              legend.direction = "horizontal") +
        scale_fill_manual(values=val_bar_fill) +
        guides(fill=g) +
        labs(fill="") + coord_cartesian(clip = "off") +
        annotate("text", label = gene_type_, 
                 x = ceiling(length(upset_data$sizes$exclusive_intersection)/2)+0.5, 
                 y = (max(upset_data$sizes$exclusive_intersection)-5), size = 6)
      
    ),
    set_sizes=(upset_set_size(geom = geom_col(width=0.6, mapping=aes(fill=present, alpha=opacity)))
               +ylab('Source Size')+ xlab("Source") + scale_fill_manual(
                 values=val_bar_fill
                 
               )),
    themes=upset_modify_themes(
      list(
        'intersections_matrix'=theme(text=element_text(size=16), 
                                     axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.length.x=unit(0,'lines')),
        'overall_sizes'=theme(axis.text.x=overall_sizes_x_text,
                              axis.title.x = overall_sizes_x_label)
      )
    ),
    name="",
    wrap=TRUE,
    height_ratio=1.5,
    stripes=upset_stripes(
      geom=geom_segment(size=5),
      colors=unname(sets_color[sets])
    ),
    add = add,
    limits = c(-900,900)
  )
  
}


p_bar <- ggplot(df_bar, aes(x=source, y=value, fill=present, pattern = remove)) +
  ggpattern::geom_bar_pattern(
    position="stack", stat="identity",
    color = "white",
    pattern_angle = 45
    ) +
  theme_minimal() +
  theme(
    axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.text.x=element_text(size=14),
    axis.title.x = element_blank(),
    strip.text = element_text(size=18),
    plot.margin = margin(-5, 0, 0, 0, "lines"),
    legend.position = "none"
  ) +
  scale_fill_manual(values=val_bar_fill) +
  ggpattern::scale_pattern_manual(values=c('none','stripe')) +
  labs(y="# of Alleles", x="") +
  facet_wrap(.~gene_type, scales="free_y", ncol=7)

p <- wrap_elements(
  Reduce(f='/', p_upset) + plot_layout(heights = c(1,5/6,5/6)) & theme(plot.margin = margin(1, 0, 0, 0, "lines"))
)

p_all <- wrap_elements(
  p / p_bar + plot_layout(heights = c(1,1/8)) + plot_annotation(tag_levels = "A") & 
    theme(plot.tag = element_text(size = 25),plot.margin = margin(0.5, 0, 0, 0, "lines"))
)


ggsave(p_all, filename = "figures/figure3.pdf", width = 14, height = 16, units = "in", dpi = 300)

for(gene_type_ in c("IGHD","IGHJ","IGLJ","IGKJ")){
    df <- df_macaque_subset[gene_type==gene_type_]
    ## remove any columns that are all values == FALSE
    df <- df[, lapply(.SD, function(col) if (all(col == FALSE)) NULL else col)]
    sets <- intersect(intersect_, colnames(df))
    upset_data <- ComplexUpset::upset_data(df, intersect = sets, sort_intersections = F)
    
    if(gene_type_ %in% c("IGHD","IGHJ")){
      for(source_ in c("IMGT","KIMDB","Trios")){
        n1 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection==sprintf("%s-Baseline", source_)& upset_data$with_sizes$in_database==FALSE & upset_data$with_sizes$exclusive_intersection==sprintf("%s-Baseline", source_),])
        n2 <- nrow(upset_data$with_sizes[upset_data$with_sizes$intersection==sprintf("%s-Baseline", source_)& upset_data$with_sizes$in_database==TRUE & upset_data$with_sizes$exclusive_intersection==sprintf("%s-Baseline", source_),])
        print(sprintf("%s %s not in database: %s", gene_type_, source_, n1))
        print(sprintf("%s %s in database: %s", gene_type_, source_, n2))
        print(sprintf("%s %s total: %s", gene_type_, source_, n2+n1))
      }
    }
    print("----------------------------------")
    
    order_ <- sort(upset_data$sizes$exclusive_intersection[!grepl("Novel",names(upset_data$sizes$exclusive_intersection))], decreasing = F)
    upset_data$sorted$intersections <- c("Novel",names(order_))
    upset_data$sorted$groups <- sets
    g <- guide_legend(title="")
    g <- if(gene_type_ != "IGHD") "none"
    overall_sizes_x_label <- if(gene_type_ %in% c("IGKV","IGKJ")) element_text(size=16)  else element_blank()
    overall_sizes_x_text <- if(gene_type_ %in% c("IGKV","IGKJ")) element_text(size=16)  else element_blank()
    add <- if(gene_type_ %in% c("IGHV","IGHD")) TRUE else FALSE
    p_upset[[gene_type_]] <- upset_adapted(
    data = upset_data,
    base_annotations=list(
    '# of Alleles'=intersection_size2(
    counts=TRUE,
    mapping=aes(fill=present),
    text_colors = c("black","black")
    )+theme(axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    legend.text = element_text(size=20),
    legend.title = element_text(size=21),
    #legend.position = c(0.1, 1.3),
    legend.position = "top",
    legend.direction = "horizontal") +
    scale_fill_manual(values=val_bar_fill, drop=FALSE) +
    guides(fill=g) +
    labs(fill="") + coord_cartesian(clip = "off") +
    annotate("text", label = gene_type_,
    x = ceiling(length(upset_data$sizes$exclusive_intersection)/2),
    y = (max(upset_data$sizes$exclusive_intersection)), size = 6)
    ),
    set_sizes=(upset_set_size(geom = geom_col(width=0.6, mapping=aes(fill=present, alpha=opacity)))
    +ylab('Source Size')+ xlab("Source") + scale_fill_manual(
    values=val_bar_fill
    )),
    themes=upset_modify_themes(
    list(
    'intersections_matrix'=theme(text=element_text(size=16),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length.x=unit(0,'lines')),
    'overall_sizes'=theme(axis.text.x=overall_sizes_x_text,
    axis.title.x = overall_sizes_x_label)
    )
    ),
    name="",
    wrap=TRUE,
    height_ratio=1.5,
    stripes=upset_stripes(
    geom=geom_segment(size=5),
    colors=unname(sets_color[sets])
    ),
    add = add,
    limits = c(-75,15)
    ) 
}
p <- wrap_elements(
Reduce(f='/', p_upset[c("IGHD","IGHJ","IGLJ","IGKJ")]) + plot_layout(ncol = 1)
)
ggsave(p, filename = "figures/figure3_supp.pdf", width = 14, height = 18, units = "in", dpi = 300)



