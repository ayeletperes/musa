## infer ASC for the MUSA set
## functions
max_ <- function(x, inf.rm =  T) {
  if(length(x) == 0){
    return(NA)
  }
  return(max(x[is.finite(x)], na.rm = T))
}

min_ <- function(x, inf.rm =  T) {
  if(length(x) == 0){
    return(NA)
  }
  return(min(x[is.finite(x)], na.rm = T))
}

which.max_ <- function(x, inf.rm =  T) {
  if(length(x) == 0){
    return(NA)
  }
  return(which.max(x[is.finite(x)]))
}

# Compute the distance matrix between sequences.
# @param sequences A character vector of sequences.
# @param method The method to compute the distance. Either "lv" for Levenshtein distance or "hamming" for Hamming distance.
# @param trimm_prim_3 The number of nucleotides to trim from the 3' end of the sequences.
# @param type The type of output. Either "both" for both the dis and a matrix, "matrix" or "dist.
compute_reference_distance <- function(sequences, method = "lv", trimm_prim_3 = NULL, type = "both", quiet = FALSE){
  
  ## asset that trimm_prim_3 is either null or numeric
  if(!is.null(trimm_prim_3) && !is.numeric(trimm_prim_3)){
    stop("trimm_prim_3 must be either NULL or numeric")
  }
  ## assert that trimm_prim_3 is larger than 1
  if(!is.null(trimm_prim_3) && trimm_prim_3 < 1){
    stop("trimm_prim_3 must be larger than 1")
  }
  ## assert that the method is either lv or hamming
  if(!(method %in% c("lv", "hamming"))){
    stop("method must be either lv or hamming")
  }
  
  ## if trimm_prim_3 is not null, trimm the sequences
  if(!is.null(trimm_prim_3)){
    sequences <- substr(sequences, 1, trimm_prim_3)
  }
  
  ## if method is hamming make sure that all sequences have equal length. if not add N to max length and through a message
  if(method == "hamming"){
    if(!quiet) message(sprintf("sequences have been padded with N to the max length (%s).", max(nchar(sequences))))
    sequences <- gsub("\\s", "N", format(sequences, width = max(nchar(sequences))))
  }
  
  ## compute the distance matrix
  dist <- stringdist::stringdistmatrix(sequences, method = method, useNames = "names")
  
  if(type == "both"){
    return(list(dist = dist, matrix = as.matrix(dist)))
  }else{
    if(type == "matrix"){
      return(as.matrix(dist))
    }else{
      return(dist)
    }
  }
}

getNClusters <- function(g, n_cluster, range_min = 0, range_max = 6, max_steps = 20, method = "louvain") {
  # Initialize binary search variables
  this_step <- 0
  this_min <- range_min
  this_max <- range_max
  pre_min_cluster <- 0
  pre_max_cluster <- NULL
  min_update <- FALSE
  max_update <- FALSE
  
  # Store the best partition and number of clusters found so far
  best_partition <- NULL
  closest_clusters <- NULL
  
  while (this_step < max_steps) {
    this_step <- this_step + 1
    
    # Determine the current resolution
    if (this_step == 1) {
      this_resolution <- this_min + (this_max - this_min) / 2
    } else {
      if (max_update && !min_update) {
        this_resolution <- this_min + (this_max - this_min) * (n_cluster - pre_min_cluster) / (pre_max_cluster - pre_min_cluster)
      } else if (min_update && !max_update) {
        if (!is.null(pre_max_cluster)) {
          this_resolution <- this_min + (this_max - this_min) * (n_cluster - pre_min_cluster) / (pre_max_cluster - pre_min_cluster)
        } else {
          this_resolution <- this_min + (this_max - this_min) / 2
        }
      }
    }
    
    # Perform clustering
    if (method == "louvain") {
      partition <- cluster_louvain(g,resolution = this_resolution, 
                                   weights = E(g)$weight)
    } else if (method == "leiden") {
      partition <- cluster_leiden(g, objective_function = "CPM",
                                  resolution = this_resolution, 
                                  weights = E(g)$weight)
    } else {
      stop("Error: Unsupported method. Choose 'louvain' or 'leiden'.")
    }
    
    # Determine the number of clusters
    this_clusters <- length(unique(partition$membership))
    
    # Update search bounds
    if (this_clusters > n_cluster) {
      this_max <- this_resolution
      pre_max_cluster <- this_clusters
      min_update <- FALSE
      max_update <- TRUE
    } else if (this_clusters < n_cluster) {
      this_min <- this_resolution
      pre_min_cluster <- this_clusters
      min_update <- TRUE
      max_update <- FALSE
    } else {
      # Found exact match
      closest_clusters <- this_clusters
      best_partition <- partition
      break
    }
    
    # Store the closest match so far
    if (is.null(closest_clusters) || abs(this_clusters - n_cluster) < abs(closest_clusters - n_cluster)) {
      closest_clusters <- this_clusters
      best_partition <- partition
    }
  }
  
  # Return the best partition found
  return(list(partition = best_partition, clusters = closest_clusters, best_resolution = this_resolution))
}

optimize_silhouette_and_get_clusters_macaque <- function(g, distance_matrix, clusters_merged_df, target_clusters = 80, 
                                                         resoltuion_range_low = 0.1, 
                                                         resolution_range_high = 0.5, 
                                                         max_steps = 20, ncors = -1) {
  if(ncors==-1){
    ncors = detectCores() - 1
  }
  # Initialize data storage for plotting
  results <- data.frame(
    Resolution = numeric(),
    ClusterCount = integer(),
    SilhouetteLeiden = numeric(),
    RandIndex = numeric(),
    Calinski_Harabasz = numeric()
  )
  
  result_leiden <- getNClusters(g, n_cluster = target_clusters, range_min = 0, range_max = 1, max_steps = max_steps, method = "leiden")
  
  min_resolution <- result_leiden$best_resolution-result_leiden$best_resolution*resoltuion_range_low
  max_resolution <- result_leiden$best_resolution*resolution_range_high+result_leiden$best_resolution
  resolution_range <- seq(min_resolution, max_resolution, by=0.05)
  if(!result_leiden$best_resolution %in% resolution_range){
    resolution_range <- c(seq(min_resolution, max_resolution, by=0.05),result_leiden$best_resolution)
  }
  
  partitions <- list()
  # Loop through each cluster count in the target range
  pb <- progress:: progress_bar$new(
    format = "(:spin) [:bar] :percent eta: :eta",
    total = length(resolution_range), clear = FALSE, width = 60)
  
  for (resolution in resolution_range) {
    
    cluster_memberships <- mclapply(1:100, function(i) {
      set.seed(i)
      partition <- cluster_leiden(g, objective_function = "CPM",resolution = resolution, weights = E(g)$weight)
      rm(.Random.seed, envir=globalenv())
      g_membership <- unlist(membership(partition))
      paste0(paste0(names(g_membership), '#', g_membership),collapse = ";")
    }, mc.cores = ncors, mc.set.seed = F, mc.cleanup = T)
    
    
    membership_counts <- table(unlist(cluster_memberships))
    most_frequent_membership <- names(which.max(membership_counts))
    most_frequent_count <- max(membership_counts)
    g_membership_optimum <- unlist(lapply(
      strsplit(unlist(strsplit(most_frequent_membership,";")),"#"),
      function(x){
        setNames(as.integer(x[2]),x[1])
      }))
    n_clusters <- length(unique(g_membership_optimum))
    partitions[[paste0("res-",resolution)]] <- g_membership_optimum
    
    df2 <- data.table(allele = names(g_membership_optimum), leiden = g_membership_optimum)
    df_merge <- merge(clusters_merged_df, df2, by = "allele")
    
    if (n_clusters > 1) {
      sil_score_leiden <- silhouette(g_membership_optimum, distance_matrix)
      if(!any(is.na(sil_score_leiden))){
        sil_score_leiden <- mean(sil_score_leiden[, 3])
      }
      ri_leiden <- rand.index(df_merge$leiden, df_merge$hc)
      ch_leiden <- intCriteria(distance_matrix, as.integer(g_membership_optimum),c("Calinski_Harabasz"))$calinski_harabasz
    } else {
      sil_score_leiden <- -1  # Set to -1 for invalid clustering
      ri_leiden <- NA
      ch_leiden <- NA
    }
    # Record scores in results data frame
    results <- rbind(results, data.frame(
      Resolution = resolution,
      ClusterCount = n_clusters,
      SilhouetteLeiden = sil_score_leiden,
      RandIndex = ri_leiden,
      Calinski_Harabasz = ch_leiden
    ))
    pb$tick()
  }
  
  results$Calinski_Harabasz_norm <- (results$Calinski_Harabasz - min_(results$Calinski_Harabasz)) / (max_(results$Calinski_Harabasz) - min_(results$Calinski_Harabasz))
  
  # Plotting
  p <- ggplot(results, aes(x = Resolution)) +
    geom_line(aes(y = SilhouetteLeiden, color = "Silhouette")) +
    geom_line(aes(y = Calinski_Harabasz_norm, color = "Calinski Harabasz")) +
    geom_line(aes(y = RandIndex, color = "Rand Index")) +
    geom_vline(xintercept = result_leiden$best_resolution, linetype = "dashed", color = "palegreen2") +
    labs(
      x = "Resolution",
      y = "Value",
      color = "Metric"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("Calinski Harabasz" = "lightblue2", "Silhouette" = "lightsalmon", 
                                  "Rand Index" = "tomato"))
  
  max_silhouette <- which.max(results$SilhouetteLeiden)
  best_resolution <- results$Resolution[max_silhouette]
  best_partition <- partitions[[max_silhouette]]
  
  # Return the best result
  return(list(
    results = results,
    partitions = partitions,
    plot = p,
    best_resolution = best_partition,
    best_partition = best_partition
  ))
}


message_log <- function(event = "run completed", project = "r-server", channel = "runs", icon = "🔥", notify = TRUE) {
  token <- readRDS("token.rds")
  
  headers = c(
    'Content-Type' = 'application/json',
    'Authorization' = paste0('Bearer ', token)
  )
  
  body <- list(
    project = project,
    channel = channel,
    event = event,
    icon = icon,
    notify = notify
  )
  
  body <- rjson::toJSON(body)
  
  res <- httr::VERB("POST", url = "https://api.logsnag.com/v1/log", body = body, add_headers(headers))
  
}

if (!require("librarian")) install.packages("librarian")

librarian::shelf(
  ComplexHeatmap,
  tigger,
  data.table,
  stringdist,
  pbapply,
  ggplot2,
  patchwork,
  ggrepel,
  dplyr,
  tidyr,
  ggbreak,
  GetoptLong,
  grid,
  igraph,
  cluster,
  mclust,
  fpc,
  parallel,
  clusterCrit,
  fossil,
  Biostrings,
  ggExtra,
  colorspace,
  pracma,
  ggtree,
  httr,
  gridtext,
  ggforce
)

### cluster the asc of the macaque.

setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

MUSA_set <- fread("data/MUSA_data/MUSA_set_all_2025-02-04_with_iglabel.csv")

day <- as.Date(format(Sys.time(), "%Y-%m-%d")) + 1 #2024-12-04 # latest stable version - "2024-12-04"
# add a day so it won't clash with the previous version

lapply(c("IGH","IGK","IGL"), function(chain) dir.create(file.path("../macaque_asc/work/",chain,day), showWarnings = F))
## copy the database to the work directory
lapply(c("IGH","IGK","IGL"), function(chain) file.copy(sprintf("data/MUSA_data/db/macaca_mulatta_%s_with_trio.csv",tolower(chain)), file.path("../macaque_asc/work/",chain,day), overwrite = T))

## create the new IgLabel data

lapply(c("IGH","IGK","IGL"), function(chain){
  lapply(c("V","J"), function(segment){
    g <- paste0(chain,segment)
    sequences <- MUSA_set[gene_type==g,.(iglabel,seq_gapped)]
    
    if(segment=="V"){
      tigger::writeFasta(setNames(sequences$seq_gapped,sequences$iglabel), 
                         file.path("../macaque_asc/work",chain,day,paste0(segment,"_gapped_iglabel.fasta")))
    }
    
    tigger::writeFasta(setNames(gsub("[.]","",sequences$seq_gapped),sequences$iglabel), 
                       file.path("../macaque_asc/work",chain,day,paste0(segment,"_iglabel.fasta")))
  })
  
  if(chain=="IGH"){
    sequences <- MUSA_set[gene_type=="IGHD",.(iglabel,seq_gapped)]
    tigger::writeFasta(setNames(gsub("[.]","",sequences$seq_gapped),sequences$iglabel), 
                       file.path("../macaque_asc/work",chain,day,"D_iglabel.fasta"))
  }
})

ref <- grep("V_iglabel[.]fasta",grep("_iglabel.fasta",grep(day,list.files("../macaque_asc/work","fasta", full.names = T, recursive = T),value=T), value=T),invert=T, value=T)
names(ref) <- paste0("macaque_",gsub("_gapped","",gsub("_iglabel","",paste0(basename(dirname(dirname(ref))),tools::file_path_sans_ext(basename(ref))))))
ref <- ref[(!grepl("IG[KL]D", names(ref)))]
ref <- sapply(ref, readIgFasta)

## compute the distance matrix for each.
ref_distance <- pbsapply(names(ref), function(segment){
  method <- ifelse(grepl("V", segment), "hamming", "lv")
  compute_reference_distance(ref[[segment]], method = method, trimm_prim_3 = NULL, type = "both")
}, USE.NAMES = T, simplify = F)

ref_cluster_complete <- pbsapply(names(ref_distance), function(segment){
  hclust(ref_distance[[segment]]$dist, method = "complete")
}, USE.NAMES = T, simplify = F)

ref_families_cluster <- pbsapply(names(ref_distance), function(segment){
  cutoff <- max(nchar(gsub("[.]","",ref[[segment]]))) * 0.25
  dendextend::cutree(ref_cluster_complete[[segment]], h = cutoff, order_clusters_as_data = T)
}, USE.NAMES = T, simplify = F)

segments <- names(ref)
ref_families_cluster_update <- list()
for(segment in segments){
  g <- gsub("macaque_","",segment)
  
  ## first for the cluster families we should retain the closest number of the seed set naming. 
  tmp <- MUSA_set[gene_type==g]
  tmp2 <- ref_families_cluster[[segment]]
  tmp[,family_cluster:=tmp2[iglabel]]
  tmp[,seed_family:=gsub(g,"",alakazam::getFamily(allele))]
  seed_families <- gsub(g,"",unique(alakazam::getFamily(tmp$allele)))
  print(segment)
  print(seed_families)
  print("----------")
  fam_number <- setNames(rep(0,length(seed_families)), as.character(sort(as.integer(seed_families))))
  # get expressed families
  expressed_clusters <- tmp[,.(present = any(present)), by = family_cluster][present==TRUE,family_cluster]
  non_expressed_clusters <- unique(tmp2)[!unique(tmp2) %in% expressed_clusters]
  
  ## first we find the closest cluster to the seed family
  for(fam in as.character(sort(as.integer(seed_families)))){
    # extracting all related families
    cls <- tmp[seed_family==fam,unique(family_cluster)]
    # keep just expressed clusters
    cls <- cls[cls %in% expressed_clusters]
    # counting the number of alleles connected to the seed family
    cls_alleles <- tapply(tmp[family_cluster %in% cls,seed_family], 
                          tmp[family_cluster %in% cls,family_cluster], function(x) table(x)[[fam]], simplify = F)
    # choosing the cluster with the most alleles in the seed
    cl_index <- names(cls_alleles)[which.max(cls_alleles)]
    # attributing the cluster index
    if(fam_number[[fam]]==0 & !is.null(cl_index)){
      if(!as.numeric(cl_index) %in% fam_number){
         fam_number[[fam]] <- as.numeric(cl_index)
      }
    }
  }
  
  ## special correction for IGHJ
  
  if(g=="IGHJ"){
    fam_number[[as.character(non_expressed_clusters)]] <- non_expressed_clusters
    fam_number_complete <- fam_number
  }else{
    ## remove any clusters above the first zero
    fam_number_complete <- fam_number[1:(min(which(fam_number==0))-1)]
    
    ## Second, we will enumerate first families that have at least one expressed alleles.
    reamining_clusters <- unique(tmp2)[!unique(tmp2) %in% fam_number_complete]
    reamining_expressed <- expressed_clusters[expressed_clusters %in% reamining_clusters]
    reamining_not_expressed <- non_expressed_clusters[non_expressed_clusters %in% reamining_clusters]
    
    ## attribute the expressed families
    next_family <- max(as.integer(names(fam_number_complete[fam_number_complete!=0]))) + 1
    for(cl in reamining_expressed){
      fam_number_complete[[as.character(next_family)]] <- cl
      next_family <- next_family + 1
    }
    
    ## attribute the non expressed
    for(cl in non_expressed_clusters){
      fam_number_complete[[as.character(next_family)]] <- cl
      next_family <- next_family + 1
    }
    
  }
  
  clust_to_fam <- setNames(names(fam_number_complete),as.character(fam_number_complete))
  
  ref_families_cluster_update_df <- data.table(imgt_allele = names(ref_families_cluster[[segment]]), family_cluster = ref_families_cluster[[segment]])
  ref_families_cluster_update_df[,Family:=clust_to_fam[as.character(family_cluster)]]
  ref_families_cluster_update[[segment]] <- ref_families_cluster_update_df
}
## check the number of alleles for each

lapply(ref, length)

load_data <- FALSE
if(load_data){
  load(paste0("MUSA/MUSA_data_",day,".rds"))
}else{
  save(ref, ref_distance, ref_cluster_complete, ref_families_cluster_update, file = paste0("MUSA/MUSA_data_",day,".rds"))  
}

## cluster the sets; this should be run in the server

server <- FALSE
if(server){
  #load(paste0("MUSA/MUSA_data_",day,".rds"))
  results_lieden_optimized_summarised_macaque <- list()
  segments <- names(ref)
  for(segment in segments){
    distance_matrix <- ref_distance[[segment]]$matrix
    complete_threshold <- ifelse(grepl("V", segment), 0.05, 0.1)
    N <- max(nchar(gsub("[.]","",ref[[segment]]))) * complete_threshold
    # normalize the distance matrix by the max
    distance_matrix_norm <- distance_matrix/max(distance_matrix)
    similarity_matrix <- -log(distance_matrix_norm)/max_(-log(distance_matrix_norm))
    similarity_matrix[!is.finite(similarity_matrix)] <- 0
    g <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "undirected", weighted = TRUE, diag = FALSE)
    hc <- hclust(as.dist(distance_matrix), method = "complete")
    complete_clusters <- cutree(hc, h = N)
    target_clusters <- length(unique(complete_clusters))
    cluster_range <- round(target_clusters*0.25)
    imgt <- setNames(as.numeric(factor(alakazam::getGene(names(ref[[segment]]), strip_d = F, omit_nl = F))), names(ref[[segment]])) 
    complete_clusters_df <- data.table(allele = names(complete_clusters), hc = unlist(complete_clusters))
    imgt_clusters_df <- data.table(allele = names(imgt), imgt = unlist(imgt))
    clusters_merged_df <- merge.data.table(complete_clusters_df, imgt_clusters_df, by = "allele")
    results_lieden_optimized_summarised_macaque[[segment]] <- optimize_silhouette_and_get_clusters_macaque(g, distance_matrix, 
                                                                                                           clusters_merged_df = clusters_merged_df, 
                                                                                                           target_clusters = target_clusters)
    
    p <- results_lieden_optimized_summarised_macaque[[segment]]$plot
    ggsave(filename = sprintf("MUSA/figures/%s_%s_seed_smooth_log_silhouette.png", day, segment), 
           plot = p, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
    
    g_membership_optimum <- results_lieden_optimized_summarised_macaque[[segment]]$best_partition
    
    l <- ifelse(N>=10,5,2)
    s <- seq(0,N+3,by=l)
    colors <- c("firebrick","gold2", "cornflowerblue","pink")
    colors <- colors[1:(length(s)-1)]
    colors <- c(colors,"white")
    col_fun = circlize::colorRamp2(s, colors)
    or <- colnames(distance_matrix)[hc$order]
    
    col_leiden <- setNames(sapply(seq_along(unique(g_membership_optimum[colnames(distance_matrix)])), function(i){
      if(i %% 2 == 0){
        return("black")
      }else{
        return("white")
      }
    }),unique(g_membership_optimum[or]))
    
    col_complete <- setNames(sapply(seq_along(unique(complete_clusters[colnames(distance_matrix)])), function(i){
      if(i %% 2 == 0){
        return("black")
      }else{
        return("white")
      }
    }), unique(complete_clusters[or]))
    
    column_ha = HeatmapAnnotation(
      leiden = factor(g_membership_optimum[colnames(distance_matrix)]),
      complete = factor(complete_clusters[colnames(distance_matrix)]),
      show_legend = c(FALSE, FALSE),
      col = list(leiden = col_leiden, complete = col_complete)
    )
    
    d <- ComplexHeatmap::Heatmap(distance_matrix, col = col_fun,
                                 cluster_rows = hc,
                                 cluster_columns = hc,
                                 top_annotation = column_ha,
                                 show_row_names = FALSE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 name = "Substitutions")
    
    pdf(sprintf("data/MUSA_data/figures/%s_%s_seed_smooth_log_shilhouette.pdf",day, segment), width = 35, height = 15)
    draw(d)
    dev.off()
    #message_log(event = sprintf("Complete %s asc inference", segment))
  }
  save(results_lieden_optimized_summarised_macaque, file = "data/MUSA_data/results_lieden_optimized_summarised_macaque.rda")
}else{
  load("data/MUSA_data/results_lieden_optimized_summarised_macaque.rda")  
}

## rename the clusters. 
work_dir <- '/home/ayelet/BIU Dropbox/Ayelet Peres/ig_allele_papers_test_and_share/macaque_asc/work'
distributed_dir <- '/home/ayelet/BIU Dropbox/Ayelet Peres/ig_allele_papers_test_and_share/macaque_asc/distributed/2025-01-03'
db_file <- function(chain, day) sprintf('/home/ayelet/BIU Dropbox/Ayelet Peres/ig_allele_papers_test_and_share/macaque_asc/work/%s/%s/macaca_mulatta_%s_with_trio.csv', chain, day, tolower(chain))
day <- as.Date(format(Sys.time(), "%Y-%m-%d")) + 1

python_path <- '/home/ayelet/anaconda3/bin/python'
segments <- names(ref)
for(segment in segments){
  
  clusters_asc <- data.table(
    imgt_allele = names(results_lieden_optimized_summarised_macaque[[segment]]$best_partition),
    Allele_Cluster = results_lieden_optimized_summarised_macaque[[segment]]$best_partition
  )
  
  clusters_family <- ref_families_cluster_update[[segment]][,.(imgt_allele, Family)]

  clusters_table <- merge.data.table(
    clusters_family, clusters_asc, by = "imgt_allele", all = T
  )
  
  clusters_table[,index:=1:.N, by=.(Family,Allele_Cluster)]
  
  prefix <- strsplit(segment, "_")[[1]][2]
  
  clusters_table[,new_allele:=paste0(
    prefix,Family,"-G",Allele_Cluster,"*",ifelse(index<10,0,""),index
  )]
  
  seg <- substr(prefix,4,4)
  chain <- substr(prefix,1,3)
  path <- file.path(work_dir,chain,day)
  dir.create(path, recursive = T, showWarnings = F)
  
  # Create the optimum clusters table. Output file: *_asc.csv table
  asc_file <- sprintf("%s_asc.csv", seg)
  data.table::fwrite(clusters_table, file.path(path,asc_file), sep = ",", quote = F, row.names = F)
  
  # Run create_iglabel_cluster_mapping.py to create the mapping. Output file: *_asc_iglabel.csv
  iglabel_file <- sprintf("%s_asc_iglabel.csv", seg)
  system(sprintf("cd '%s' && %s ../../../python/create_iglabel_cluster_mapping.py %s %s", path, python_path, asc_file, iglabel_file))
  
  # Run rename_fasta.py to rename the fasta. Input files: *_asc_iglabel.csv IG*_label.fasta  Output file: IG*.fasta
  iglabel_fasta_file <- sprintf("%s_iglabel.fasta", seg)
  fasta_file <- sprintf("%s.fasta", seg)
  tigger::writeFasta(gsub("[.]","",ref[[segment]]), file.path(path, iglabel_fasta_file))
  
  system(sprintf("cd '%s' && %s ../../../python/rename_fasta.py %s %s %s", path, python_path, iglabel_file, iglabel_fasta_file, fasta_file))
  
  if(seg=="V"){
    iglabel_gapped_fasta_file <- sprintf("%s_gapped_iglabel.fasta", seg)
    fasta_gapped_file <- sprintf("%s_gapped.fasta", seg)
    tigger::writeFasta(ref[[segment]], file.path(path, iglabel_gapped_fasta_file))
    system(sprintf("cd '%s' && %s ../../../python/rename_fasta.py %s %s %s", path, python_path, iglabel_file, iglabel_gapped_fasta_file, fasta_gapped_file))
  }
  
  # Run merge_into_new_clusters.py to retain the names. Input files: old distribution *.fasta , new work directory IG*_label.fasta *_asc.csv V.fasta V_merge_mapping.csv V_merge_stats.txt , macaca_mulatta_db.csv
  
  ## to retain the new family assignment from the new majority rule, we will change the family assignment in the old distribution to the new one. keeping the asc names
  ## for the merging. 
  distributed_file <- file.path(distributed_dir,chain,sprintf("%s.fasta", seg))
  mapping_file <- sprintf("%s_merge_mapping.csv", seg)
  stats_file <- sprintf("%s_merge_stats.txt", seg)
  db_file_new_allele <- db_file(chain, day) #sprintf("macaca_mulatta_%s.csv", chain)
  
  system(sprintf("cd '%s' && %s ../../../python/merge_into_new_clusters_families_flag.py %s '%s' %s %s '%s' %s %s %s --new_families_only --fill_gaps", path, python_path, prefix, 
                 distributed_file, iglabel_fasta_file, asc_file, db_file_new_allele, fasta_file, mapping_file, stats_file))
  
}







# create the wheel for each segment

clusters_tables <- list()

segments <- names(ref)
for(segment in segments){
  prefix <- strsplit(segment, "_")[[1]][2]
  seg <- substr(prefix,4,4)
  chain <- substr(prefix,1,3)
  path <- file.path(work_dir,chain,day)
  ## get the new asc iglabel table.
  iglabel_file <- sprintf("%s_asc_iglabel.csv", seg)
  df <- fread(file.path(path, iglabel_file))
  ## get the new merged mapping
  mapping_file <- sprintf("%s_merge_mapping.csv", seg)
  df_new_merged <- fread(file.path(path, mapping_file))
  
  df_new_merged[,original_allele:=paste0(prefix,"0-",label,"*00")]
  df_new_merged <- merge.data.table(df_new_merged, df[,.(original_allele, Family, Allele_Cluster)], by = "original_allele", all.y = T)
  present_asc <- MUSA_set[gene_type==prefix & present==TRUE,iglabel]
  df_new_merged[,present:=original_allele %in% present_asc]
  df_new_merged[,new_Family:=gsub(prefix,"",alakazam::getFamily(new_name, strip_d = F, omit_nl = F))]
  ## families that all their ASCs are not expressed change to 0
  df_new_merged[,tag:=ifelse(any(present),unique(new_Family),0),by=.(new_Family)]
  df_new_merged[,new_tag:=apply(df_new_merged[,.(new_name,tag)],1,function(x){
    base <- strsplit(x[["new_name"]],"-")[[1]]
    allele <- base[2]
    paste0(prefix,x[["tag"]],"-",allele, collapse = "")
  })]
  df_new_merged[,new_tag:=gsub(" ","",new_tag)]
  clusters_tables[[segment]] <- df_new_merged
}

clusters_tables_df <- rbindlist(clusters_tables)

## save the bag of alleles with the new asc names. 
iglabe_new_allele <- setNames(clusters_tables_df$new_name, clusters_tables_df$original_allele)
iglabe_new_tag <- setNames(clusters_tables_df$new_tag, clusters_tables_df$original_allele)

MUSA_set[,new_allele:=iglabe_new_allele[iglabel]]
MUSA_set[,new_tag:=iglabe_new_tag[iglabel]]

fwrite(MUSA_set, sprintf("data/MUSA_data/MUSA_set_all_%s_with_asc.csv",day))

## load the original seed set with the novel alleles from the genomic data and create the rss reference.
seed_set <- fread("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_set_and_novel_2025-02-04_update.csv", sep = ",")

old_allele_to_new_allele <- setNames(MUSA_set$new_allele, MUSA_set$allele)
seed_set[,new_allele:=old_allele_to_new_allele[allele_update]]

old_allele_to_iglabel <- setNames(MUSA_set$iglabel, MUSA_set$allele)
seed_set[,iglabel:=old_allele_to_iglabel[allele_update]]

fwrite(seed_set, sprintf("data/rhesus_macaque_data/seed_set_with_novel_alleles/seed_set_with_novel_%s_with_asc.csv",day))

## merge the table with the MUSA set
tmp <- MUSA_set[,.(allele,chain,gene_type,iglabel,new_allele,seq,seq_gapped,
             kimdb, `rhgldb+`, imgt, trios, vrc, guo, novel, notes, 
             novel_notes, new_tag)]

columns_to_process <- c("kimdb", "rhgldb+", "imgt", "trios", "vrc", "guo", "novel", "notes", "novel_notes", "new_tag")

seed_set_all <- seed_set
for (col_name in columns_to_process) {
  new_allele_map <- setNames(MUSA_set[[col_name]], MUSA_set$new_allele)
  seed_set_all[[col_name]] <- new_allele_map[seed_set$new_allele]
}

## bind seed and MUSA

seed_alleles <- seed_set_all$new_allele
merged_data <- rbind(seed_set_all, MUSA_set[!new_allele %in% seed_alleles,.(allele,chain,gene_type,iglabel,new_allele,seq,seq_gapped,
                                                                             kimdb, `rhgldb+`, imgt, trios, vrc, guo, novel, notes, 
                                                                             novel_notes, new_tag)], fill = T)

fwrite(merged_data, sprintf("MUSA/MUSA_with_rss_and_leader_and_asc_%s.csv",day))

## filter the bag of alleles just to sure rss, give the new family and label for the rss alignment.

seed_set_rss_filter <- seed_set[!grepl("NONAMER not found|RSS not found", digger_notes, ignore.case = TRUE) & 
                                              sample_count_genomic>0 & 
                                              sample_count_AIRRseq>0 & 
                                              sure_subject_count>0,]

fwrite(seed_set_rss_filter, sprintf("data/MUSA_data/seed_set_with_novel_rss_filtered_%s_with_asc.csv",day))

## save the data needed for figure 4
save(clusters_tables, ref_distance, ref, file = "data/figure4/figure4_data.rda")
