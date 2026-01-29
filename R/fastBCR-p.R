#' Function: Preprocess datasets and check required columns of paired data
#' @description To infer clonal families successfully, the dataset should include essential columns:
#' 1. v_call_heavy: heavy chain V gene with or without allele,
#' 2. j_call_heavy: heavy chain J gene with or without allele,
#' 3. junction_aa_heavy: amino acid translation of the heavy chain junction,
#' 4. v_call_light: light chain V gene with or without allele,
#' 5. j_call_light: light chain J gene with or without allele,
#' 6. junction_aa_light: amino acid translation of the light chain junction
#' @param paired_raw_data_list A list of paired data.
#'
#' @return Processed data as input for clonal family inference
#' @export
paired.preprocess <- function(paired_raw_data_list,
                              productive_only = FALSE,
                              count_col_name = NA,
                              count_filter_thre = NA) {
  for (i in seq_along(paired_raw_data_list)) {
    # check heavy/ light chains
    if (!"v_call_heavy" %in% colnames(paired_raw_data_list[[i]]) ||
        !"j_call_heavy" %in% colnames(paired_raw_data_list[[i]]) ||
        !"junction_aa_heavy" %in% colnames(paired_raw_data_list[[i]]) ||
        !"v_call_light" %in% colnames(paired_raw_data_list[[i]]) ||
        !"j_call_light" %in% colnames(paired_raw_data_list[[i]]) ||
        !"junction_aa_light" %in% colnames(paired_raw_data_list[[i]])) {
      stop("lack essential columns(BCR heavy/light chain)!")
    }

    paired_raw_data_list[[i]]$v_call <- paired_raw_data_list[[i]]$v_call_heavy
    paired_raw_data_list[[i]]$j_call <- paired_raw_data_list[[i]]$j_call_heavy
    paired_raw_data_list[[i]]$junction_aa <- paired_raw_data_list[[i]]$junction_aa_heavy
  }
  paired_data_list <- data.preprocess(raw_data_list = paired_raw_data_list, productive_only, count_col_name, count_filter_thre)
  return(paired_data_list)
}

BCR.clusters.p<- function(input,
                          min_depth_thre = 3,
                          max_depth_thre = 1000,
                          overlap_thre = 0.1,
                          consensus_thre = 0.8) {

  input$v_call <- input$v_call_heavy
  input$j_call <- input$j_call_heavy
  input$junction_aa <- input$junction_aa_heavy

  bcr_clusters <- c()

  ### 1. Fast k-mer pre-clustering
  ## 1.1 VJ partition
  VJ <- paste(input$v_call, input$j_call, sep = "_")
  VJ_table <- sort(table(VJ), decreasing = T)
  VJ_sati <- names(VJ_table)[which(as.numeric(VJ_table) >= min_depth_thre)]
  for (vj in 1:length(VJ_sati)) {
    tmp.vj <- VJ_sati[vj]
    vj.loc <- which(VJ == tmp.vj)
    vj.seqs <- input[vj.loc, ]
    # Length pre-clustering
    L <- nchar(vj.seqs$junction_aa)
    L_table <- sort(table(L), decreasing = T)
    L_sati <- names(L_table)[which(as.numeric(L_table) >= min_depth_thre)]
    if (length(L_sati) == 0) next

    ## 1.2 k-mer clustering
    pre_clusters <- pre_clustering(L, L_sati, vj.seqs, min_depth_thre)

    ## 1.3 Merge clusters with same seed
    if (length(pre_clusters) == 0) next
    filt_clusters <- merge_clusters(pre_clusters)

    ### 2.Optimized clustering
    ## 2.1 Further combine pre-clusters
    sort_clusters <- Sort_clu(filt_clusters)
    clu.indexs <- lapply(sort_clusters, function(x) x$clonotype_index)
    df <- data.frame(
      index = 1:length(clu.indexs),
      indexs = I(clu.indexs)
    )
    df <- dplyr::distinct(df, indexs, .keep_all = T)
    clu.index <- df$index
    sort_clusters <- sort_clusters[clu.index]

    cc <- length(sort_clusters)
    if (cc == 1) {
      bcr_clusters <- c(bcr_clusters, sort_clusters)
      next
    }

    bcr_clusters <- combine_pre_clusters(bcr_clusters, sort_clusters, cc, overlap_thre)
  }

  if (length(bcr_clusters) != 0) {
    final_clusters <- list()
    for (i in seq_along(bcr_clusters)) {
      cluster <- bcr_clusters[[i]]
      light_VJ <- paste(cluster$v_call_light, cluster$j_call_light, sep = "_")
      cluster$light_VJ <- light_VJ
      split_clusters <- split(cluster, light_VJ)
      split_clusters <- lapply(split_clusters, function(df) df[, !(names(df) %in% c("light_VJ"))])
      final_clusters <- c(final_clusters, split_clusters)
    }

    # Remove temporary columns and sort clusters
    final_clusters <- Sort_clu(final_clusters)
    bcr_clusters <- final_clusters[sapply(final_clusters, nrow) >= 3]

    # 2.2 Filter candidates with low consensus score
    bcr_clusters <- Sort_clu(bcr_clusters)
    bcr_clusters_depth <- as.numeric(unlist(sapply(bcr_clusters, function(x) nrow(x))))
    depth.filt.loc <- which(bcr_clusters_depth > max_depth_thre)
    if (length(depth.filt.loc) != 0) {
      bcr_clusters <- bcr_clusters[-depth.filt.loc]
    }
    bcr_clusters_consensus <- consensus_scores(bcr_clusters)
    consensus.filt.loc <- which(bcr_clusters_consensus < consensus_thre)
    if (length(consensus.filt.loc) != 0) {
      bcr_clusters <- bcr_clusters[-consensus.filt.loc]
    }
    if(length(bcr_clusters) != 0){
      for (i in 1:length(bcr_clusters)) {
        tmp.bcr_clusters <- bcr_clusters[[i]]
        tmp.bcr_clusters <- tmp.bcr_clusters[, -which(colnames(tmp.bcr_clusters) %in% c("mode", "loc", "len", "kmer"))]
        bcr_clusters[[i]] <- tmp.bcr_clusters
      }
    }

  }

  return(bcr_clusters)
}

BCR.clusters.unfilter.p<- function(input,
                                   min_depth_thre = 3,
                                   max_depth_thre = 1000,
                                   overlap_thre = 0.1) {
  input$v_call <- input$v_call_heavy
  input$j_call <- input$j_call_heavy
  input$junction_aa <- input$junction_aa_heavy

  bcr_clusters <- c()
  min_depth_thre <- min_depth_thre + floor(nrow(input) / 1e+05)

  ### 1. Fast k-mer pre-clustering
  ## 1.1 VJ partition
  VJ <- paste(input$v_call, input$j_call, sep = "_")
  VJ_table <- sort(table(VJ), decreasing = T)
  VJ_sati <- names(VJ_table)[which(as.numeric(VJ_table) >= min_depth_thre)]
  for (vj in 1:length(VJ_sati)) {
    tmp.vj <- VJ_sati[vj]
    vj.loc <- which(VJ == tmp.vj)
    vj.seqs <- input[vj.loc, ]
    # Length pre-clustering
    L <- nchar(vj.seqs$junction_aa)
    L_table <- sort(table(L), decreasing = T)
    L_sati <- names(L_table)[which(as.numeric(L_table) >= min_depth_thre)]
    if (length(L_sati) == 0) next

    ## 1.2 k-mer clustering
    pre_clusters <- pre_clustering(L, L_sati, vj.seqs, min_depth_thre)

    ## 1.3 Merge clusters with same seed
    if (length(pre_clusters) == 0) next
    filt_clusters <- merge_clusters(pre_clusters)

    ### 2.Optimized clustering
    ## 2.1 Further combine pre-clusters
    sort_clusters <- Sort_clu(filt_clusters)
    clu.indexs <- lapply(sort_clusters, function(x) x$clonotype_index)
    df <- data.frame(
      index = 1:length(clu.indexs),
      indexs = I(clu.indexs)
    )
    df <- dplyr::distinct(df, indexs, .keep_all = T)
    clu.index <- df$index
    sort_clusters <- sort_clusters[clu.index]

    cc <- length(sort_clusters)
    if (cc == 1) {
      bcr_clusters <- c(bcr_clusters, sort_clusters)
      next
    }

    bcr_clusters <- combine_pre_clusters(bcr_clusters, sort_clusters, cc, overlap_thre)
  }

  if (length(bcr_clusters) != 0) {
    final_clusters <- list()
    for (i in seq_along(bcr_clusters)) {
      cluster <- bcr_clusters[[i]]
      light_VJ <- paste(cluster$v_call_light, cluster$j_call_light, sep = "_")
      cluster$light_VJ <- light_VJ
      split_clusters <- split(cluster, light_VJ)
      split_clusters <- lapply(split_clusters, function(df) df[, !(names(df) %in% c("light_VJ"))])
      final_clusters <- c(final_clusters, split_clusters)
    }

    # Remove temporary columns and sort clusters
    final_clusters <- Sort_clu(final_clusters)
    bcr_clusters <- final_clusters[sapply(final_clusters, nrow) >= 3]
  }

  return(bcr_clusters)
}

#' Predict Public Antibody Scores
#'
#' This function predicts the public antibody scores for input BCR sequences using pre-trained Python models.
#'
#' @param data A data frame containing BCR data. Must include columns: `cdr1`, `cdr2`, `cdr3`, and `vgene`.
#' @param model A character string specifying the prediction model to use. Options are:
#'   - `"cdrh"`: Model trained on heavy chain CDR regions.
#'   - `"cdrl"`: Model trained on light chain CDR regions.
#'   - `"cdrh3"`: Model trained on heavy chain CDR3 region only.
#'   - `"cdrl3"`: Model trained on light chain CDR3 region only.
#' @param python_env A character string specifying the Python environment to use. Defaults to `"public"`.
#'
#' @return A data frame with an additional column, `public_score`, containing the predicted public antibody scores.
#' @export
#' @import reticulate
predict_public_antibody <- function(data, model = "cdrh", python_env = "public") {
  use_condaenv(python_env, required = TRUE)

  py$data <- data
  py$model <- model

  python_code <- "
import os
import pandas as pd
import torch
from PubBCRPredictor import PubBCRPredictor_Runner, MLP
from BCR_V_BERT import BCR_V_BERT_Runner

def vgene_process(vgene):
    if '*' in vgene:
        vgene = vgene.split('*')[0]
    if '/OR' in vgene:
        vgene = vgene.replace('/OR','')
    if 'S' in vgene:
        vgene = vgene.replace('S','-')
    if vgene.count('-')>1:
        vgene = vgene.split('-')[0]+'-'+vgene.split('-')[1]
    return vgene

data['cdr'] = data['cdr1'] + '|' + data['cdr2'] + '|' + data['cdr3']
sequence = data['cdr'].values
data['vgene']= data['vgene'].apply(vgene_process)
vgenes = data['vgene'].values
cdr3s = data['cdr3'].values

model_name = model
if model == 'cdrh':
    bcr_v_bert = BCR_V_BERT_Runner(model='cdrh')
    pub_runner = PubBCRPredictor_Runner(model='cdrh')
    feature = bcr_v_bert.embed(sequence,vgenes)

elif model == 'cdrl':
    bcr_v_bert = BCR_V_BERT_Runner(model='cdrl')
    pub_runner = PubBCRPredictor_Runner(model='cdrl')
    feature = bcr_v_bert.embed(sequence,vgenes)

elif model == 'cdrh3':
    bcr_v_bert = BCR_V_BERT_Runner(model='cdrh3')
    pub_runner = PubBCRPredictor_Runner(model='cdrh3')
    feature = bcr_v_bert.embed(cdr3s, vgenes)

elif model == 'cdrl3':
    bcr_v_bert = BCR_V_BERT_Runner(model='cdrl3')
    pub_runner = PubBCRPredictor_Runner(model='cdrl3')
    feature = bcr_v_bert.embed(cdr3s, vgenes)

prob = pub_runner.predict(feature)

# 如果 prob 是 GPU Tensor, 转为 numpy 格式
if isinstance(prob, torch.Tensor):
    prob = prob.cpu().detach().numpy()
"
  reticulate::py_run_string(python_code)
  prob_in_r <- py$prob

  data$public_score <- prob_in_r
  return(data)
}

#' Filter Public Antibodies by Threshold
#'
#' This function filters predicted public antibody results into high and low public antibody categories based on a threshold.
#'
#' @param prediction_results A data frame containing prediction results, including a `public_score` column.
#' @param threshold A numeric value specifying the public score threshold for classification. Defaults to 0.5.
#'
#' @return A list containing two data frames:
#'   - `high_public`: Rows where `public_score` is greater than or equal to the threshold.
#'   - `low_public`: Rows where `public_score` is less than the threshold.
#' @export
#' @importFrom dplyr select filter
filter_public_antibodies <- function(prediction_results, threshold = 0.5) {
  high_public <- prediction_results[prediction_results$public_score >= threshold, ]
  low_public <- prediction_results[prediction_results$public_score < threshold, ]
  list(
    high_public = high_public,
    low_public = low_public
  )
}

#' zscore helper
zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s  <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mu) / s
}

#' Convert v_identity to percent if it looks like fraction
.as_percent_identity <- function(v) {
  v_num <- suppressWarnings(as.numeric(v))
  if (all(is.na(v_num))) return(v_num)
  mx <- suppressWarnings(max(v_num, na.rm = TRUE))
  # heuristic: if values are in [0,1] treat as fraction
  if (is.finite(mx) && mx <= 1.5) v_num <- v_num * 100
  v_num
}

#' Prepare input for predict_public_antibody (heavy)
.make_heavy_input <- function(df_cluster) {
  req <- c("cdr1_heavy","cdr2_heavy","cdr3_heavy","v_call_heavy")
  miss <- setdiff(req, colnames(df_cluster))
  if (length(miss) > 0) {
    stop("Heavy chain columns missing: ", paste(miss, collapse = ", "))
  }
  dplyr::transmute(
    df_cluster,
    cdr1 = .data$cdr1_heavy,
    cdr2 = .data$cdr2_heavy,
    cdr3 = .data$cdr3_heavy,
    vgene = .data$v_call_heavy
  )
}

#' Prepare input for predict_public_antibody (light)
.make_light_input <- function(df_cluster) {
  req <- c("cdr1_light","cdr2_light","cdr3_light","v_call_light")
  miss <- setdiff(req, colnames(df_cluster))
  if (length(miss) > 0) {
    stop("Light chain columns missing: ", paste(miss, collapse = ", "))
  }
  dplyr::transmute(
    df_cluster,
    cdr1 = .data$cdr1_light,
    cdr2 = .data$cdr2_light,
    cdr3 = .data$cdr3_light,
    vgene = .data$v_call_light
  )
}

#' Compute per-seq SHM from v_identity
.add_shm_cols <- function(df_cluster) {
  if (!("v_identity_heavy" %in% colnames(df_cluster))) {
    stop("Missing column: v_identity_heavy")
  }
  if (!("v_identity_light" %in% colnames(df_cluster))) {
    stop("Missing column: v_identity_light")
  }
  vh <- .as_percent_identity(df_cluster$v_identity_heavy)
  vl <- .as_percent_identity(df_cluster$v_identity_light)
  dplyr::mutate(
    df_cluster,
    shm_heavy = 100 - vh,
    shm_light = 100 - vl
  )
}

#' Main: annotate clusters with public predictions + SHM means + downstream flags
#'
#' @param cluster_list list of data.frame (fastBCR clustering output)
#' @param python_env conda env name used by predict_public_antibody (default "r-py-env")
#' @param P_cut z-score cutoff for public_score_z (default 0.82)
#' @param SHM_cut shm_mean cutoff (default 1.88)
#' @param heavy_model public model for heavy (default "cdrh")
#' @param light_model public model for light (default "cdrl")
#'
#' @return list(
#'   cluster_list = augmented cluster_list (adds public_heavy/public_light/shm_heavy/shm_light),
#'   cluster_summary = data.frame per cluster with means + zscores,
#'   df_flag = cluster_summary with is_public/public_origin labels
#' )
annotate_public_and_flag <- function(
  cluster_list,
  python_env = "r-py-env",
  P_cut = 0.82,
  SHM_cut = 1.88,
  heavy_model = "cdrh",
  light_model = "cdrl"
) {
  if (!is.list(cluster_list) || length(cluster_list) == 0) {
    stop("cluster_list must be a non-empty list of data.frames.")
  }

  # ensure needed pkgs (keep lightweight; user can attach themselves too)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Please install purrr.")
  if (!exists("predict_public_antibody", mode = "function")) {
    stop("predict_public_antibody() not found. Please load fastBCR (or the module providing it).")
  }

  cluster_ids <- names(cluster_list)
  if (is.null(cluster_ids) || any(cluster_ids == "")) {
    cluster_ids <- paste0("cluster_", seq_along(cluster_list))
  }

  # 1) per cluster: predict public heavy/light + compute shm per seq
  cluster_list_aug <- purrr::map2(cluster_list, cluster_ids, function(df_cluster, cid) {
    if (!is.data.frame(df_cluster)) {
      stop("Each element in cluster_list must be a data.frame. Problem at: ", cid)
    }

    # SHM columns
    df_cluster <- .add_shm_cols(df_cluster)

    # heavy public prediction
    heavy_in <- .make_heavy_input(df_cluster)
    pred_h <- predict_public_antibody(heavy_in, model = heavy_model, python_env = python_env)
    if (!("public_score" %in% colnames(pred_h))) {
      stop("predict_public_antibody heavy output must contain column: public_score")
    }

    # light public prediction
    light_in <- .make_light_input(df_cluster)
    pred_l <- predict_public_antibody(light_in, model = light_model, python_env = python_env)
    if (!("public_score" %in% colnames(pred_l))) {
      stop("predict_public_antibody light output must contain column: public_score")
    }

    # attach to original cluster df (row order assumed consistent)
    dplyr::mutate(
      df_cluster,
      public_heavy = as.numeric(pred_h$public_score),
      public_light = as.numeric(pred_l$public_score)
    )
  })

  # 2) per cluster: summarize means
  cluster_summary <- purrr::map2_dfr(cluster_list_aug, cluster_ids, function(df_cluster, cid) {
    dplyr::tibble(
      cluster_id = cid,
      n_seq = nrow(df_cluster),
      public_heavy_mean = mean(df_cluster$public_heavy, na.rm = TRUE),
      public_light_mean = mean(df_cluster$public_light, na.rm = TRUE),
      shm_heavy_mean = mean(df_cluster$shm_heavy, na.rm = TRUE),
      shm_light_mean = mean(df_cluster$shm_light, na.rm = TRUE)
    )
  })

  # 3) zscore + downstream flags
  df <- dplyr::mutate(
    cluster_summary,
    public_heavy_z = zscore(public_heavy_mean),
    public_light_z = zscore(public_light_mean),
    public_score_z = (public_heavy_z + public_light_z) / 2,
    shm_mean = (shm_heavy_mean + shm_light_mean) / 2
  )

  df_flag <- dplyr::mutate(
    df,
    is_public = public_score_z >= P_cut,
    public_origin = dplyr::case_when(
      is_public & shm_mean <  SHM_cut ~ "Naive-derived (filtered)",
      is_public & shm_mean >= SHM_cut ~ "Memory-derived (kept)",
      TRUE                            ~ "Non-public"
    )
  )

  list(
    cluster_list = cluster_list_aug,
    cluster_summary = df,
    df_flag = df_flag
  )
}
