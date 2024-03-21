Sort_clu <- function(clu) {
  len.lis <- sapply(clu, function(x) length((x$clonotype_index)))
  index <- order(len.lis, decreasing = T)
  sort <- vector("list", length = length(index))
  sort[c(seq_along(index))] <- clu[index]
  
  return(sort)
}

sort_msa <- function(seqs, seqs0) {
  l <- length(seqs)
  loc <- c()
  all <- c(1:l)
  all <- as.character(all)
  for (i in 1:l) {
    tmp.loc <- which(seqs[i] == seqs0)
    t <- length(tmp.loc)
    if (t > 1) {
      for (tt in 1:t) {
        tmp.loc0 <- tmp.loc[tt]
        cha <- as.character(tmp.loc0)
        if (length(grep(paste("^", cha, "$", sep = ""), all)) != 0) {
          all <- gsub(paste("^", cha, "$", sep = ""), "", all)
          loc <- c(loc, tmp.loc0)
          break
        } else {
          next
        }
      }
    } else {
      cha <- as.character(tmp.loc)
      all <- gsub(paste("^", cha, "$", sep = ""), "", all)
      loc <- c(loc, tmp.loc)
    }
  }
  
  return(loc)
}

#' Function: Calculate consensus scores of clonal families infered by fastBCR
#'
#' @param bcr_clusters Clonal families infered by fastBCR
#'
#' @return A list where each element is the consensus score of each bcr_cluster
#' @export
consensus_scores <- function(bcr_clusters) {
  con.lis <- c()
  for (i in seq_along(bcr_clusters)) {
    seqs_aa0 <- bcr_clusters[[i]]$junction_aa
    sink("default matrix.txt")
    MSAalign_aa <- msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW")
    sink()
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    splt <- strsplit(seqs_aa, "")
    each_AA <- NULL
    num <- length(splt)
    AA_len <- length(splt[[1]])
    most_AA <- c()
    for (j in 1:AA_len) {
      AA_lis <- sapply(splt, function(x) x[j])
      each_AA[[j]] <- AA_lis
      most_AA <- c(most_AA, names(sort(table(AA_lis), decreasing = TRUE))[1])
    }
    gap.loc <- which(most_AA == "-")
    if (length(gap.loc) == 0) {
      con <- sapply(each_AA, function(x) as.numeric(sort(table(x), decreasing = T))[1] / num)
    } else {
      con <- c()
      for (k in seq_along(each_AA)) {
        if (any(gap.loc == k)) {
          tmp.con <- as.numeric(sort(table(each_AA[[k]]), decreasing = T))[2] / num
          con <- c(con, tmp.con)
        } else {
          tmp.con <- as.numeric(sort(table(each_AA[[k]]), decreasing = T))[1] / num
          con <- c(con, tmp.con)
        }
      }
    }
    con.mean <- mean(con)
    con.lis <- c(con.lis, con.mean)
  }
  
  return(con.lis)
}

pre_clustering <- function(L, L_sati, vj.seqs, cluster_thre){
  pre_clusters <- NULL
  pc <- 1
  # first three k-mers
  start1 <- 3
  end1 <- 7
  start2 <- 4
  end2 <- 8
  start3 <- 5
  end3 <- 9
  for (l in 1:length(L_sati)) {
    tmp.l <- L_sati[l]
    tmp.seqs <- vj.seqs[which(L == tmp.l), ]
    tmp.junc <- tmp.seqs$junction_aa
    # last three k-mers
    half <- ceiling(as.numeric(tmp.l) / 2)
    start4 <- half - 1
    end4 <- half + 3
    start5 <- half
    end5 <- half + 4
    start6 <- half + 1
    end6 <- half + 5
    
    # k-mer indexing
    k_mer1 <- substr(tmp.junc, start1, end1)
    k_mer2 <- substr(tmp.junc, start2, end2)
    k_mer3 <- substr(tmp.junc, start3, end3)
    k_mer4 <- substr(tmp.junc, start4, end4)
    k_mer5 <- substr(tmp.junc, start5, end5)
    k_mer6 <- substr(tmp.junc, start6, end6)
    
    # IGHJ6 optimization
    j_call <- unique(tmp.seqs$j_call)
    if (j_call == "IGHJ6") {
      j6.kmer <- c("YYYYG", "YYYGM", "YYGMD", "YGMDV", "GMDVW")
      kmer4.none <- unique(c(which(k_mer4 %in% j6.kmer), grep("..YYY", k_mer4), grep(".YYYG", k_mer4)))
      kmer5.none <- unique(c(which(k_mer5 %in% j6.kmer), grep("..YYY", k_mer5), grep(".YYYG", k_mer5)))
      kmer6.none <- unique(c(which(k_mer6 %in% j6.kmer), grep("..YYY", k_mer6), grep(".YYYG", k_mer6)))
      k_mer4[kmer4.none] <- ""
      k_mer5[kmer5.none] <- ""
      k_mer6[kmer6.none] <- ""
    }
    
    k_mer.lis <- c(k_mer1, k_mer2, k_mer3, k_mer4, k_mer5, k_mer6)
    index <- rep(tmp.seqs$clonotype_index, times = 6)
    k_mer.fre <- sort(table(k_mer.lis), decreasing = T)
    kk <- which(as.numeric(k_mer.fre) >= cluster_thre)
    k_mer.clu <- names(k_mer.fre)[kk]
    empty.kmer <- which(k_mer.clu == "")
    if (length(empty.kmer) != 0) {
      k_mer.clu <- k_mer.clu[-empty.kmer]
    }
    
    for (mm in k_mer.clu) {
      index.loc <- which(k_mer.lis == mm)
      indexs <- index[index.loc]
      indexs <- unique(indexs)
      tmp.loc <- which(tmp.seqs$clonotype_index %in% indexs)
      if (length(tmp.loc) < cluster_thre) {
        next
      }
      tmp.dd <- tmp.seqs[tmp.loc, ]
      tmp.dd$kmer <- mm
      
      # filtering k-mer position
      cdr3s <- tmp.dd$junction_aa
      kmer <- tmp.dd$kmer[1]
      kmer.start <- stringr::str_locate(cdr3s, kmer)[, 1]
      tmp.table <- table(kmer.start)
      tmp.max <- max(tmp.table)
      kmer.start.mode <- mean(as.numeric(names(which(tmp.table == tmp.max))))
      
      left <- which(abs(kmer.start - kmer.start.mode) < 1.5)
      if (length(left) < cluster_thre) {
        next
      }
      kmer.start <- kmer.start[left]
      tmp.table <- table(kmer.start)
      tmp.max <- max(tmp.table)
      kmer.start.mode <- mean(as.numeric(names(which(tmp.table == tmp.max))))
      tmp.dd$mode <- kmer.start.mode
      tmp.dd$loc <- pc
      tmp.dd$len <- unique(nchar(cdr3s))
      pre_clusters[pc] <- list(tmp.dd[left, ])
      pc <- pc + 1
    }
  }
  
  return(pre_clusters)
}

merge_clusters <- function(pre_clusters){
  all_kmer <- sapply(pre_clusters, function(x) unique(x$kmer))
  kmer_table <- table(all_kmer)
  potential_kmer <- names(kmer_table)[which(as.numeric(kmer_table) > 1)]
  pk <- length(potential_kmer)
  
  if (pk == 0) {
    filt_clusters <- pre_clusters
  } else {
    filt_clusters <- pre_clusters[-which(all_kmer %in% potential_kmer)]
    for (i in 1:pk) {
      tmp.kmer <- potential_kmer[i]
      kmer.loc <- which(all_kmer == tmp.kmer)
      kmer.dd <- pre_clusters[kmer.loc]
      sati.length <- sapply(kmer.dd, function(x) unique(x$len))
      for (j in 1:length(kmer.loc)) {
        clu.length <- unique(kmer.dd[[j]]$len)
        clu.loc <- unique(kmer.dd[[j]]$loc)
        clu.mode <- unique(kmer.dd[[j]]$mode)
        
        plus1 <- clu.length + 1
        minus1 <- clu.length - 1
        plus1.grep <- length(grep(paste("^", plus1, "$", sep = ""), sati.length))
        minus1.grep <- length(grep(paste("^", minus1, "$", sep = ""), sati.length))
        
        if (plus1.grep > 0 && minus1.grep == 0) {
          clu.plus.loc <- unique(kmer.dd[[which(sati.length == plus1)]]$loc)
          clu.plus.mode <- unique(kmer.dd[[which(sati.length == plus1)]]$mode)
          if (abs(clu.mode - clu.plus.mode) < 1) {
            clu.loc <- c(clu.loc, clu.plus.loc)
          }
        } else if (plus1.grep == 0 && minus1.grep > 0) {
          clu.minus.loc <- unique(kmer.dd[[which(sati.length == minus1)]]$loc)
          clu.minus.mode <- unique(kmer.dd[[which(sati.length == minus1)]]$mode)
          if (abs(clu.mode - clu.minus.mode) < 1) {
            clu.loc <- c(clu.loc, clu.minus.loc)
          }
        } else if (plus1.grep > 0 && minus1.grep > 0) {
          clu.plus.loc <- unique(kmer.dd[[which(sati.length == plus1)]]$loc)
          clu.plus.mode <- unique(kmer.dd[[which(sati.length == plus1)]]$mode)
          clu.minus.loc <- unique(kmer.dd[[which(sati.length == minus1)]]$loc)
          clu.minus.mode <- unique(kmer.dd[[which(sati.length == minus1)]]$mode)
          if (abs(clu.mode - clu.plus.mode) < 1) {
            clu.loc <- c(clu.loc, clu.plus.loc)
          }
          if (abs(clu.mode - clu.minus.mode) < 1) {
            clu.loc <- c(clu.loc, clu.minus.loc)
          }
        } else {
          
        }
        if (length(clu.loc) > 1) {
          filt.clu <- pre_clusters[[clu.loc[1]]]
          for (l in 2:length(clu.loc)) {
            filt.clu <- rbind(filt.clu, pre_clusters[[clu.loc[l]]])
          }
          filt.clu <- filt.clu[!duplicated(filt.clu$clonotype_index), ]
        } else {
          filt.clu <- pre_clusters[[clu.loc]]
        }
        filt.len <- mean(nchar(filt.clu$junction_aa))
        filt.clu$len <- filt.len
        filt_clusters <- c(filt_clusters, list(filt.clu))
      }
    }
  }
  
  return(filt_clusters)
}

combine_pre_clusters <- function(bcr_clusters, sort_clusters, cc, overlap_thre){
  clu_info <- vector(mode = "list", length = cc)
  for (i in 1:cc) {
    clu_info[[i]] <- list(
      Id = sort_clusters[[i]]$clonotype_index,
      len = unique(sort_clusters[[i]]$len),
      index = i
    )
  }
  clu.index <- c(1:cc)
  dd <- 1
  while (T) {
    # Candidate
    clu1.id <- clu_info[[dd]]$Id
    clu1.len <- clu_info[[dd]]$len
    
    dd.else <- dd + 1
    clu.else.len <- sapply(clu_info, function(x) x$len)[dd.else:cc]
    clu.else.index <- unlist(sapply(clu_info, function(x) x$index)[dd.else:cc])
    satis.len <- which(abs(clu.else.len - clu1.len) <= 1)
    satis.index <- clu.else.index[satis.len]
    
    comb.loc <- c()
    if (length(satis.index) != 0) {
      for (inx in satis.index) {
        tmp.loc <- which(clu.index == inx)
        clu2.id <- clu_info[[tmp.loc]]$Id
        clu2.size <- length(clu2.id)
        per <- length(intersect(clu1.id, clu2.id)) / clu2.size
        if (per >= overlap_thre) {
          comb.loc <- c(comb.loc, tmp.loc)
        }
      }
    }
    
    n_num <- length(comb.loc)
    if (n_num != 0) {
      dd.id <- clu_info[[dd]]$Id
      dd.index <- clu_info[[dd]]$index
      dd.len <- clu_info[[dd]]$len
      for (n in comb.loc) {
        tmp.id <- clu_info[[n]]$Id
        tmp.index <- clu_info[[n]]$index
        tmp.len <- clu_info[[n]]$len
        dd.id <- c(dd.id, tmp.id)
        dd.index <- c(dd.index, tmp.index)
        dd.len <- c(dd.len, tmp.len)
      }
      dd.id <- unique(dd.id)
      clu_info[[dd]]$Id <- dd.id
      clu_info[[dd]]$index <- dd.index
      clu_info[[dd]]$len <- mean(dd.len)
      clu_info <- clu_info[-comb.loc]
      clu.index <- clu.index[-comb.loc]
    } else {
      dd.inx <- clu_info[[dd]]$index
      com.clu <- sort_clusters[[dd.inx[1]]]
      com.clu <- as.data.frame(com.clu)
      dd.inx <- dd.inx[-1]
      if (length(dd.inx) != 0) {
        for (i in dd.inx) {
          tmp <- sort_clusters[[i]]
          tmp <- as.data.frame(tmp)
          com.clu <- rbind(com.clu, tmp)
        }
        com.clu <- com.clu[!duplicated(com.clu$clonotype_index), ]
      }
      bcr_clusters <- c(bcr_clusters, list(com.clu))
      dd <- dd + 1
    }
    cc <- length(clu_info)
    if (dd == cc) {
      dd.inx <- clu_info[[dd]]$index
      com.clu <- sort_clusters[[dd.inx[1]]]
      com.clu <- as.data.frame(com.clu)
      dd.inx <- dd.inx[-1]
      if (length(dd.inx) != 0) {
        for (i in dd.inx) {
          tmp <- sort_clusters[[i]]
          tmp <- as.data.frame(tmp)
          com.clu <- rbind(com.clu, tmp)
        }
        com.clu <- com.clu[!duplicated(com.clu$clonotype_index), ]
      }
      bcr_clusters <- c(bcr_clusters, list(com.clu))
      break
    }
  }
  
  return(bcr_clusters)
}

#' Function: The core step of fastBCR to infer clonal families from processed data.
#'
#' BCR clonal family inference consists of two steps, fast k-mer-based pre-clustering and optimized clustering, to infer clonal families fast and accurately.
#'
#' @param input Processed data
#' @param cluster_thre Minimal clustering criteria. Defaluts to 3. For high efficiency, the threshold is increased by 1 for every 100,000 entries of input data.
#' @param overlap_thre The overlap coefficient threshold for merging two clusters whose selectable range is (0,1). Defaluts to 0.1. Lower thresholds may lead to overclustering while higher may lead to split of clonal families.
#' @param consensus_thre The consensus score threshold for filtering candidates. Defaluts to 0.8. A higher threshold means stricter inference of the cluster.
#'
#' @return Clonal families inferred by fastBCR
#' @export
BCR.cluster <- function(input, cluster_thre = 3,
                        overlap_thre = 0.1,
                        consensus_thre = 0.8) {
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
  
  return(bcr_clusters)
}

#' Function: Fast clonal family inference from preprocessed data
#'
#' @param pro_data_list A list where each element is the preprocessed data named after its filename
#' @param cluster_thre Minimal clustering criteria. Defaults to 3. For high efficiency, the threshold is increased by 1 for every 100,000 entries of input data.
#' @param overlap_thre The overlap coefficient threshold for merging two clusters, selectable range (0,1). Defaults to 0.1. Lower thresholds may lead to overclustering while higher thresholds may lead to the split of clonal families.
#' @param consensus_thre The consensus score threshold for filtering candidates. Defaults to 0.8. A higher threshold means stricter inference of the cluster.
#'
#' @return A list where each element is a list of clonal families inferred by fastBCR from a single sample
#' @export
data.BCR.clusters <- function(pro_data_list, cluster_thre = 3, overlap_thre = 0.1, consensus_thre = 0.8) {
  BCR_clusters_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(pro_data_list), style = 3)
  for (i in seq_along(pro_data_list)) {
    data <- pro_data_list[[i]]
    var_name <- names(pro_data_list)[i]
    processed_data <- BCR.cluster(data, cluster_thre, overlap_thre, consensus_thre)
    BCR_clusters_list[[var_name]] <- processed_data
    setTxtProgressBar(pb, i)
  }
  
  cat(paste0(
    paste0("\nYou have clustered ", length(BCR_clusters_list), " samples\n"),
    paste(paste0("Sample '", names(BCR_clusters_list), "' contains ", sapply(BCR_clusters_list, function(x) length(x)), " clonal families"), collapse = "\n")
  ))
  
  
  return(BCR_clusters_list)
}

#' Function: The core step of fastBCR to infer clonal families without consensus filtering from processed data.
#'
#' @param input Processed data
#' @param cluster_thre Minimal clustering criteria. Defaluts to 3. For high efficiency, the threshold is increased by 1 for every 100,000 entries of input data.
#' @param overlap_thre The overlap coefficient threshold for merging two clusters whose selectable range is (0,1). Defaluts to 0.1. Lower thresholds may lead to overclustering while higher may lead to split of clonal families.
#'
#' @return Unfiltered clonal families inferred by fastBCR
#' @export
BCR.cluster.unfilter <- function(input, cluster_thre = 3,
                                 overlap_thre = 0.1){
  bcr_clusters <- c()
  cluster_thre <- cluster_thre + floor(nrow(input) / 1e+05)
  
  # VJ partition
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
    
    ### Fast k-mer pre-clustering
    pre_clusters <- pre_clustering(L, L_sati, vj.seqs, cluster_thre)
    
    # Merge clusters with same seed
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
  bcr_clusters <- Sort_clu(bcr_clusters)
  
  for (i in 1:length(bcr_clusters)) {
    tmp.bcr_clusters <- bcr_clusters[[i]]
    tmp.bcr_clusters <- tmp.bcr_clusters[, -which(colnames(tmp.bcr_clusters) %in% c("mode", "loc", "len", "kmer"))]
    bcr_clusters[[i]] <- tmp.bcr_clusters
  }
  
  return(bcr_clusters)
}
