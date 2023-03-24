Sort_clu <- function(clu) {
  len.lis <- sapply(clu, function(x) length((x$sequence_id)))
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

similar_con <- function(bcr_clusters) {
  con.lis <- c()
  for (i in seq_along(bcr_clusters)) {
    seqs_aa0 <- bcr_clusters[[i]]$junction_aa
    MSAalign_aa <- msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW")
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

BuildBCRlineage <- function(input, cluster_thre) {
  junction_aa <- input$junction_aa
  len_aa <- nchar(junction_aa)
  half <- ceiling(len_aa / 2)

  # 5-mer*6
  start1 <- 3
  end1 <- 7
  start2 <- 4
  end2 <- 8
  start3 <- 5
  end3 <- 9
  start4 <- half - 1
  end4 <- half + 3
  start5 <- half
  end5 <- half + 4
  start6 <- half + 1
  end6 <- half + 5

  k_mer1 <- substr(junction_aa, start1, end1)
  k_mer2 <- substr(junction_aa, start2, end2)
  k_mer3 <- substr(junction_aa, start3, end3)
  k_mer4 <- substr(junction_aa, start4, end4)
  k_mer5 <- substr(junction_aa, start5, end5)
  k_mer6 <- substr(junction_aa, start6, end6)

  # IGHJ6 optimization
  j_call <- input$j_call
  j6.loc <- which(j_call == "IGHJ6")
  if (length(j6.loc) != 0) {
    j6.kmer <- c("YYYYG", "YYYGM", "YYGMD", "YGMDV", "GMDVW")
    j6.kmer4.loc <- unique(c(which(k_mer4 %in% j6.kmer), grep("..YYY", k_mer4), grep(".YYYG", k_mer4)))
    kmer4.none <- intersect(j6.loc, j6.kmer4.loc)
    j6.kmer5.loc <- unique(c(which(k_mer5 %in% j6.kmer), grep("..YYY", k_mer5), grep(".YYYG", k_mer5)))
    kmer5.none <- intersect(j6.loc, j6.kmer5.loc)
    j6.kmer6.loc <- unique(c(which(k_mer6 %in% j6.kmer), grep("..YYY", k_mer6), grep(".YYYG", k_mer6)))
    kmer6.none <- intersect(j6.loc, j6.kmer6.loc)
    k_mer4[kmer4.none] <- ""
    k_mer5[kmer5.none] <- ""
    k_mer6[kmer6.none] <- ""
  }

  k_mer.lis <- c(k_mer1, k_mer2, k_mer3, k_mer4, k_mer5, k_mer6)
  id <- rep(input$sequence_id, times = 6)

  k_mer.fre <- sort(table(k_mer.lis), decreasing = T)
  kk <- which(as.numeric(k_mer.fre) >= cluster_thre)
  k_mer.clu <- names(k_mer.fre)[kk]
  empty.kmer <- which(k_mer.clu == "")
  if (length(empty.kmer) != 0) {
    k_mer.clu <- k_mer.clu[-empty.kmer]
  }
  BCRlineage <- vector("list", length = length(k_mer.clu))
  BCRl <- 1

  for (mm in k_mer.clu) {
    id.loc <- which(k_mer.lis == mm)
    ids <- id[id.loc]
    ids <- unique(ids)
    tmp.loc <- which(input$sequence_id %in% ids)
    if (length(tmp.loc) < cluster_thre) {
      next
    }
    tmp.dd0 <- input[tmp.loc, ]

    vjcall <- paste(tmp.dd0$v_call, tmp.dd0$j_call, sep = "_")
    vjtable <- table(vjcall)
    vj <- which(as.numeric(vjtable) >= cluster_thre)
    if (length(vj) == 0) {
      next
    }
    for (v in vj) {
      tmp.vjcall <- names(vjtable)[v]
      tmp.dd <- tmp.dd0[vjcall == tmp.vjcall, ]
      tmp.dd$kmer <- mm
      BCRlineage[BCRl] <- list(tmp.dd)
      BCRl <- BCRl + 1
    }
  }

  lap <- lapply(BCRlineage, function(x) length(x$sequence_id))
  filt <- which(lap == 0)
  BCRlineage <- BCRlineage[-filt]
  return(BCRlineage)
}

#' Function: The core step of fastBCR to infer clonal families from processed data.
#'
#' BCR Clustering consists of two steps, k-mer-based pre-clustering and optimized clustering, to infer clonal families fast and accurately.
#'
#' @param input AIRR format data after Data_Processing
#' @param similarity_thre The similarity threshold for merging two clusters whose selectable range is (0,1]. Defaluts to 0.1. Lower thresholds may lead to overclustering while higher may lead to split of clonal families.
#' @param cluster_thre Minimal clustering criteria. Defaluts to (3 + floor(nrow(input)/100000).
#' @param compactness_thre The compactness threshold for filtering clusters. Defaluts to 0.8. A higher threshold means stricter inference of the cluster.
#'
#' @return Clonal families inferred by fastBCR
#' @export
#'
#' @examples
#' data("input")
#' bcr_clusters <- BCR.cluster(input)
BCR.cluster <- function(input, similarity_thre = 0.1,
                        cluster_thre = 3 + floor(nrow(input) / 100000),
                        compactness_thre = 0.8) {
  # k-mer pre-clustering
  pre_clusters <- BuildBCRlineage(input, cluster_thre)
  pre.n <- length(pre_clusters)

  # filtering:length & k-mer position
  clu.size <- sapply(pre_clusters, function(x) nchar(x$junction_aa))
  clu.table <- sapply(pre_clusters, function(x) table(nchar(x$junction_aa)))

  filt_clusters <- vector(mode = "list", length = pre.n)
  fi <- 1
  for (t in 1:pre.n) {
    sati <- which(as.numeric(clu.table[[t]]) >= cluster_thre)
    sati.length <- as.numeric(names(clu.table[[t]])[sati])
    ss <- length(sati.length)
    if (ss == 0) {
      next
    }

    left.clu <- c()
    for (ff in 1:ss) {
      clu.length <- sati.length[ff]
      filt.loc <- which(clu.size[[t]] == sati.length[ff])
      cdr3s <- pre_clusters[[t]]$junction_aa[filt.loc]
      kmer <- pre_clusters[[t]]$kmer[1]

      tmp.table <- table(stringr::str_locate(cdr3s, kmer)[, 1])
      tmp.max <- max(tmp.table)
      kmer.start.mode <- mean(as.numeric(names(which(tmp.table == tmp.max))))
      kmer.start <- stringr::str_locate(cdr3s, kmer)[, 1]
      left <- which(abs(kmer.start - kmer.start.mode) < 1.5)
      filt.loc <- filt.loc[left]
      if (length(filt.loc) < cluster_thre) {
        next
      }
      cdr3s <- cdr3s[left]
      tmp.table <- table(stringr::str_locate(cdr3s, kmer)[, 1])
      tmp.max <- max(tmp.table)
      kmer.start.mode <- mean(as.numeric(names(which(tmp.table == tmp.max))))
      len_loc <- list(
        loc = filt.loc,
        len = clu.length,
        mode = kmer.start.mode
      )
      left.clu <- c(left.clu, list(len_loc))
    }

    lc <- length(left.clu)
    if (lc == 0) {
      next
    }
    sati.length <- sapply(left.clu, function(x) x$len)

    for (f in 1:lc) {
      clu.length <- left.clu[[f]]$len
      clu.loc <- left.clu[[f]]$loc
      clu.mode <- left.clu[[f]]$mode

      plus1 <- clu.length + 1
      minus1 <- clu.length - 1
      plus1.grep <- length(grep(paste("^", plus1, "$", sep = ""), sati.length))
      minus1.grep <- length(grep(paste("^", minus1, "$", sep = ""), sati.length))

      if (plus1.grep > 0 && minus1.grep == 0) {
        clu.plus.loc <- left.clu[[which(sati.length == plus1)]]$loc
        clu.plus.mode <- left.clu[[which(sati.length == plus1)]]$mode
        if (abs(clu.mode - clu.plus.mode) < 1) {
          clu.loc <- c(clu.loc, clu.plus.loc)
        }
      } else if (plus1.grep == 0 && minus1.grep > 0) {
        clu.minus.loc <- left.clu[[which(sati.length == minus1)]]$loc
        clu.minus.mode <- left.clu[[which(sati.length == minus1)]]$mode
        if (abs(clu.mode - clu.minus.mode) < 1) {
          clu.loc <- c(clu.loc, clu.minus.loc)
        }
      } else if (plus1.grep > 0 && minus1.grep > 0) {
        clu.plus.loc <- left.clu[[which(sati.length == plus1)]]$loc
        clu.plus.mode <- left.clu[[which(sati.length == plus1)]]$mode
        clu.minus.loc <- left.clu[[which(sati.length == minus1)]]$loc
        clu.minus.mode <- left.clu[[which(sati.length == minus1)]]$mode

        if (abs(clu.mode - clu.plus.mode) < 1) {
          clu.loc <- c(clu.loc, clu.plus.loc)
        }
        if (abs(clu.mode - clu.minus.mode) < 1) {
          clu.loc <- c(clu.loc, clu.minus.loc)
        }
      } else {

      }
      filt.clu <- pre_clusters[[t]][clu.loc, ]
      filt.jun <- filt.clu$junction_aa
      filt.len <- mean(nchar(filt.jun))
      filt.clu$len <- filt.len
      filt <- list(filt.clu)
      filt_clusters[fi] <- filt
      fi <- fi + 1
    }
  }
  filt_clusters <- filt_clusters[-which(lapply(filt_clusters, function(x) length(x$sequence_id)) == 0)]

  # Sort clusters & unique
  sort_clusters <- Sort_clu(filt_clusters)
  clu.ids <- sapply(sort_clusters, function(x) x$sequence_id)
  df <- data.frame(
    index = seq_along(clu.ids),
    ids = I(clu.ids)
  )
  df <- dplyr::distinct(df, ids, .keep_all = T)
  clu.index <- df$index
  sort_clusters <- sort_clusters[clu.index]

  # Dynamic programming
  cc <- length(sort_clusters)
  clu_info <- vector(mode = "list", length = cc)
  for (i in 1:cc) {
    clu_info[[i]] <- list(
      Id = sort_clusters[[i]]$sequence_id,
      len = unique(sort_clusters[[i]]$len),
      index = i
    )
  }
  clu.index <- c(1:cc)
  bcr_clusters <- vector(mode = "list", length = cc)
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
        if (per >= similarity_thre) {
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
        com.clu <- com.clu[!duplicated(com.clu$sequence_id), ]
      }
      bcr_clusters[[dd]] <- com.clu
      dd <- dd + 1
    }
    cc <- length(clu_info)
    if (dd >= cc) {
      break
    }
  }
  bcr_clusters <- bcr_clusters[-which(lapply(bcr_clusters, function(x) length(x$sequence_id)) == 0)]
  bcr_clusters <- Sort_clu(bcr_clusters)
  con <- similar_con(bcr_clusters)
  filt.loc <- which(con < compactness_thre)
  if (length(filt.loc) != 0) {
    bcr_clusters <- bcr_clusters[-filt.loc]
  }
  return(bcr_clusters)
}
