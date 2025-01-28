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
paired.preprocess <- function(paired_raw_data_list) {
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
  paired_data_list <- data.preprocess(data_list = paired_raw_data_list)
  return(paired_data_list)
}

BCR.clusters.p<- function(input, cluster_thre = 3,
                          overlap_thre = 0.1,
                          consensus_thre = 0.8) {

  input$v_call <- input$v_call_heavy
  input$j_call <- input$j_call_heavy
  input$junction_aa <- input$junction_aa_heavy

  bcr_clusters <- c()
  cluster_thre <- cluster_thre + floor(nrow(input) / 1e+05)

  ### 1. Fast k-mer pre-clustering
  ## 1.1 VJ partition
  VJ <- paste(input$v_call, input$j_call, sep = "_")
  VJ_table <- sort(table(VJ), decreasing = T)
  VJ_sati <- names(VJ_table)[which(as.numeric(VJ_table) >= cluster_thre)]
  for (vj in 1:length(VJ_sati)) {
    tmp.vj <- VJ_sati[vj]
    vj.loc <- which(VJ == tmp.vj)
    vj.seqs <- input[vj.loc, ]
    # Length pre-clustering
    L <- nchar(vj.seqs$junction_aa)
    L_table <- sort(table(L), decreasing = T)
    L_sati <- names(L_table)[which(as.numeric(L_table) >= cluster_thre)]
    if (length(L_sati) == 0) next

    ## 1.2 k-mer clustering
    pre_clusters <- pre_clustering(L, L_sati, vj.seqs, cluster_thre)

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
    con <- consensus_scores(bcr_clusters)
    filt.loc <- which(con < consensus_thre)
    if (length(filt.loc) != 0) {
      bcr_clusters <- bcr_clusters[-filt.loc]
    }

    for (i in 1:length(bcr_clusters)) {
      tmp.bcr_clusters <- bcr_clusters[[i]]
      tmp.bcr_clusters <- tmp.bcr_clusters[, -which(colnames(tmp.bcr_clusters) %in% c("mode", "loc", "len", "kmer"))]
      bcr_clusters[[i]] <- tmp.bcr_clusters
    }

  }

  return(bcr_clusters)
}

BCR.clusters.unfilter.p<- function(input, cluster_thre = 3,
                                                 overlap_thre = 0.1) {
  input$v_call <- input$v_call_heavy
  input$j_call <- input$j_call_heavy
  input$junction_aa <- input$junction_aa_heavy

  bcr_clusters <- c()
  cluster_thre <- cluster_thre + floor(nrow(input) / 1e+05)

  ### 1. Fast k-mer pre-clustering
  ## 1.1 VJ partition
  VJ <- paste(input$v_call, input$j_call, sep = "_")
  VJ_table <- sort(table(VJ), decreasing = T)
  VJ_sati <- names(VJ_table)[which(as.numeric(VJ_table) >= cluster_thre)]
  for (vj in 1:length(VJ_sati)) {
    tmp.vj <- VJ_sati[vj]
    vj.loc <- which(VJ == tmp.vj)
    vj.seqs <- input[vj.loc, ]
    # Length pre-clustering
    L <- nchar(vj.seqs$junction_aa)
    L_table <- sort(table(L), decreasing = T)
    L_sati <- names(L_table)[which(as.numeric(L_table) >= cluster_thre)]
    if (length(L_sati) == 0) next

    ## 1.2 k-mer clustering
    pre_clusters <- pre_clustering(L, L_sati, vj.seqs, cluster_thre)

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

data.BCR.clusters.p <- function(pro_data_list, cluster_thre = 3,
                                overlap_thre = 0.1, consensus_thre = 0.8) {
  BCR_clusters_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(pro_data_list), style = 3)
  for (i in seq_along(pro_data_list)) {
    #  is_in_n <- i %in% noclust_idx
    #  if(is_in_n)
    #    next
    data <- pro_data_list[[i]]
    var_name <- names(pro_data_list)[i]
    cat("\nclustering:",i,var_name)
    processed_data <- BCR.clusters.p(data, cluster_thre, overlap_thre, consensus_thre)
    if(is.null(processed_data)){
      cat("\n", i, var_name, " has no clusterg")
      next
    }

    BCR_clusters_list[[var_name]] <- processed_data
    setTxtProgressBar(pb, i)
  }

  cat(paste0(
    paste0("\nYou have clustered ", length(BCR_clusters_list), " samples\n"),
    paste(paste0("Sample '", names(BCR_clusters_list), "' contains ", sapply(BCR_clusters_list, function(x) length(x)), " clonal families"), collapse = "\n")
  ))

  return(BCR_clusters_list)
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

data = r.data
data['cdr'] = data['cdr1'] + '|' + data['cdr2'] + '|' + data['cdr3']
sequence = data['cdr'].values
vgenes = data['vgene'].values
cdr3s = data['cdr3'].values

model_name = model
if model == 'cdrh':
    bcr_v_bert = BCR_V_BERT_Runner(model='cdrh_old')
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
