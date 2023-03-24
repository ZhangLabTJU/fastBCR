#' Function: Pieplot of V/J usage frequency
#'
#' Basic statistics of repertoire gene segment usage.
#'
#' @param data BCR repertoire
#' @param colname "v_call" for V gene or "j_call" for J gene
#'
#' @return pie chart (The top ten most frequent genes are shown, and the rest are represented by 'others')
#' @export
#'
#' @examples
#' data("input")
#' pie.freq(input, "v_call")
pie.freq <- function(data, colname) {
  gene_table <- as.data.frame(table(data[, colname]))
  gene_table <- gene_table[order(gene_table$Freq, decreasing = TRUE), ]
  gene_table$Var1 <- as.character(gene_table$Var1)
  pie.number <- gene_table[, 2]
  pie.pct <- round(pie.number / sum(pie.number) * 100)
  end <- 10
  if (length(pie.pct) > end) {
    other <- gene_table[(end + 1):nrow(gene_table), ]
    other_num <- sum(other[, 2])
    others <- c("Others", other_num)
    tmp_table <- gene_table[1:end, ]
    tmp_table <- rbind(tmp_table, others)
  } else {
    end <- length(pie.pct)
    tmp_table <- gene_table[1:end, ]
  }
  pie.number <- tmp_table$Freq
  pie.number <- as.numeric(pie.number)
  pie.pct <- round(pie.number / sum(pie.number) * 100)
  pie.pct <- paste("(", pie.pct, "%)", sep = "")
  tmp_table$Var1 <- paste(tmp_table$Var1, pie.pct, collapse = NULL)
  tmp_table$Freq <- as.numeric(tmp_table$Freq)

  pie <- ggplot(tmp_table, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", title = colname, fill = "") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(lineheight = .8, size = 45, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(face = "bold", size = 26),
      axis.title = element_text(face = "bold", size = 18, colour = "black"),
      axis.text = element_text(face = "bold", size = 16, colour = "black"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_fill_manual(breaks = tmp_table$Var1, values = RColorBrewer::brewer.pal(length(pie.number), "Paired")) +
    guides(fill = guide_legend(byrow = TRUE))
  return(pie)
}

#' Function: Histogram of junction length
#'
#' Basic statistics of repertoire junction length distribution.
#' @param data BCR repertoire
#' @param colname "DNA" for junction length or "AA" for junction_aa length
#'
#' @return junction (DNA/AA) length distribution
#' @export
#'
#' @examples
#' data("input")
#' junc.len(input, "AA")
junc.len <- function(data, colname) {
  if (colname == "DNA") {
    nchar_junc <- nchar(data[, "junction"])
  } else if (colname == "AA") nchar_junc <- nchar(data[, "junction_aa"])
  junc <- as.data.frame(nchar_junc)

  ggplot(junc, aes(x = nchar_junc)) +
    geom_histogram(
      binwidth = 2, colour = "black",
      aes(y = after_stat(density), fill = after_stat(count))
    ) +
    stat_function(fun = dnorm, color = "darkblue", size = 2, args = list(mean = mean(nchar_junc), sd = sd(nchar_junc))) +
    labs(x = "Length", y = "Density") +
    ggtitle(paste("Junction ", colname, " Length Distribution", sep = "")) +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      plot.title = element_text(lineheight = .8, size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold", size = 18, colour = "black"),
      axis.text = element_text(face = "bold", size = 15, colour = "black"),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 15)
    )
}

#######################################################################
name <- function(bcr_clusters, n, loc) {
  id <- bcr_clusters[[n]][["sequence_id"]]
  id <- id[loc]
  c_call <- bcr_clusters[[n]][["c_call"]]
  if (all(c(is.null(c_call), all(is.na(c_call))))) {
    c_call <- c_call[loc]
    c_call <- Iso.fst(c_call)
  }
  nn <- length(id)
  names <- c()
  for (i in 1:nn) {
    if (all(c(is.null(c_call), all(is.na(c_call))))) {
      tmp.c <- c_call[i]
    } else {
      tmp.c <- ""
    }
    if (i < 10) {
      tmp.name <- paste(paste("0", i, sep = ""), tmp.c, collapse = NULL)
    } else {
      tmp.name <- paste(i, tmp.c, collapse = NULL)
    }
    names <- c(names, tmp.name)
  }
  return(names)
}

Iso.fst <- function(iso) {
  for (is in seq_along(iso)) {
    is.spl <- as.vector(unlist(strsplit(iso[is], ",")))
    len <- length(is.spl)
    if (len == 1) {
      next
    }
    iso[is] <- is.spl[1]
  }
  return(iso)
}

sort.msa <- function(seqs, seqs0) {
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

#' Function: Plot of evolutionary tree and multiple sequences alignment (MSA)
#'
#' Perform MSA on all sequences in the clonal family and plot an evolutionary tree based on Hamming distance.
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param n select a cluster to plot the evolutionary tree
#'
#' @return evolutionary tree & MSA
#' @export
#'
#' @examples
#' data("bcr_clusters")
#' msa.tree(bcr_clusters, 20)
msa.tree <- function(bcr_clusters, n) {
  seqs0 <- bcr_clusters[[n]][["junction"]]
  MSAalign <- msa(Biostrings::DNAStringSet(seqs0), "ClustalW")
  seqs <- as.character(attributes(MSAalign)$unmasked)
  SEQs <- Biostrings::DNAStringSet(seqs)
  seqs_no <- gsub("-", "", seqs)
  loc <- sort_msa(seqs_no, seqs0)
  SEQs@ranges@NAMES <- name(bcr_clusters, n, loc)
  d <- as.dist(Biostrings::stringDist(SEQs, method = "hamming") / width(SEQs)[1])
  tree <- ape::bionj(d)
  p <- ggtree::ggtree(tree, options(ignore.negative.edge = TRUE)) +
    ggtree::geom_tiplab()
  d <- tidy_msa(SEQs)
  p +
    ggtree::geom_facet(geom = geom_msa, data = d, panel = "MSA", color = "Clustal") +
    ggtree::xlim_tree(0.1)
}

#' Function: Visualization of junction_aa sequences within a cluster
#'
#' Perform MSA on the junction amino acid of all sequences in the clonal family for visualization.
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param n select a cluster to visualize
#'
#' @return cluster sequences visualization.From top to bottom are title(Vgene_Jgene_junction_aaLength), seqlogo, ggmsa, and msabar
#' @export
#'
#' @examples
#' data("bcr_clusters")
#' msa.logo(bcr_clusters, 20)
msa.logo <- function(bcr_clusters, n) {
  seqs_aa0 <- bcr_clusters[[n]][["junction_aa"]]
  MSAalign_aa <- msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW")
  seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
  SEQs_aa <- Biostrings::AAStringSet(seqs_aa)
  seqs_aa_no <- gsub("-", "", seqs_aa)
  loc <- sort_msa(seqs_aa_no, seqs_aa0)
  SEQs_aa@ranges@NAMES <- name(bcr_clusters, n, loc)
  v_call <- unique(bcr_clusters[[n]][["v_call"]])
  j_call <- unique(bcr_clusters[[n]][["j_call"]])
  len <- unique(nchar(seqs_aa))
  vjlen <- paste(paste(v_call, collapse = ","), "_", paste(j_call, collapse = ","), "_Length:", len, sep = "")
  clu.title <- vjlen
  ggmsa(SEQs_aa, char_width = 0.6, seq_name = T) +
    ggtitle(label = clu.title) +
    theme(
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0),
      axis.title = element_text(size = 12, colour = "black"),
      axis.text = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(hjust = 0)
    ) +
    geom_seqlogo(adaptive = FALSE) +
    geom_msaBar()
}

SeqDist <- function(x, y) {
  if (nchar(x) != nchar(y)) {
    return(NA)
  }
  x.l <- unlist(strsplit(x, ""))
  y.l <- unlist(strsplit(y, ""))
  nn1 <- length(which(x.l != y.l))
  nn2 <- length(which(x.l == "-" & y.l != "-"))
  nn3 <- length(which(x.l != "-" & y.l == "-"))
  return(nn1 - nn2 - nn3)
}

#' Function: Calculation of SHM ratios in different isotypes in each sample
#'
#' Calculate SHM ratios in different isotypes of all clusters inferred from a repertoire.
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return four SHM ratio in different isotypes ("IGHD", "IGHM", "IGHA", and "IGHG") calculated by all clusters inferred from a repertoire
#' @export
#'
#' @examples
#' data("bcr_clusters")
#' SHM.iso(bcr_clusters)
SHM.iso <- function(bcr_clusters) {
  isotypes <- c("IGHD", "IGHM", "IGHA", "IGHG")
  all.dist <- rep(0, 4)
  all.len <- rep(0, 4)
  for (kk in bcr_clusters) {
    iso <- substr(kk[["c_call"]], 1, 4)
    seqs0 <- kk[["junction"]]
    for (is in 1:4) {
      tmp.loc <- which(iso == isotypes[is])
      if (length(tmp.loc) < 2) {
        next
      }
      tmp.seqs0 <- seqs0[tmp.loc]
      MSAalign <- msa(Biostrings::DNAStringSet(tmp.seqs0), "ClustalW")
      seqs <- as.character(attributes(MSAalign)$unmasked)
      msa.l <- unique(nchar(seqs))
      nn <- length(seqs)
      clu.l <- msa.l * nn
      dist <- matrix(0, nn, nn)
      for (ii in 1:nn) {
        for (jj in ii:nn) {
          if (jj == ii) {
            next
          }
          tmp.dist <- SeqDist(seqs[ii], seqs[jj])
          dist[jj, ii] <- dist[jj, ii] + tmp.dist
        }
      }
      dist[dist > 1] <- 0
      dist.is <- sum(dist)
      all.dist[is] <- all.dist[is] + dist.is
      all.len[is] <- all.len[is] + clu.l
    }
  }
  all.ratio <- all.dist / all.len
  return(all.ratio)
}

#' Function: Boxplot of SHM ratios in four isotypes in different samples
#'
#' Plot boxplot of SHM ratios in four isotypes in different samples. Statistical significances are evaluated
#' @param df A data.frame formed by binding SHM ratios in four isotypes  (calculated by 'SHM.iso(bcr_clusters)') from multiple samples
#'
#' @return A boxplot showing SHM ratios in four isotypes ("IGHD", "IGHM", "IGHA", "IGHG"). Statistical significances are evaluated with the Wilcoxon rank-sum test (two-sided `*P <= 0.05; **P <= 0.01; ***P <= 0.001; ****P <= 0.001`). For brevity, ns (non-significant) is not shown.
#' @export
#'
#' @examples
#' data("SHM_ratio_1")
#' data("SHM_ratio_2")
#' data("SHM_ratio_3")
#' data("SHM_ratio_4")
#' data("SHM_ratio_5")
#' df <- data.frame(
#'   Isotypes = factor(c("IGHD", "IGHM", "IGHA", "IGHG"),
#'     levels = c("IGHD", "IGHM", "IGHA", "IGHG")
#'   ),
#'   Sample1 = SHM_ratio_1,
#'   Sample2 = SHM_ratio_2,
#'   Sample3 = SHM_ratio_3,
#'   Sample4 = SHM_ratio_4,
#'   Sample5 = SHM_ratio_5
#' )
#' SHM.sample(df)
SHM.sample <- function(df) {
  df <- reshape2::melt(df)
  res <- ggpubr::compare_means(formula = value ~ Isotypes, df) %>% dplyr::filter(p.signif != "ns")
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in seq_len(nrow(res))) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = df, aes(
    x = factor(Isotypes, levels = c("IGHD", "IGHM", "IGHA", "IGHG")), y = value,
    fill = factor(Isotypes, levels = c("IGHD", "IGHM", "IGHA", "IGHG"))
  )) +
    stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6), color = "black") +
    geom_boxplot(aes(fill = factor(Isotypes, levels = c("IGHD", "IGHM", "IGHA", "IGHG"))), position = position_dodge(0.6), linewidth = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = factor(Isotypes, levels = c("IGHD", "IGHM", "IGHA", "IGHG"))), position = position_jitterdodge(), shape = 21, size = 1, alpha = 0.9) +
    labs(x = "", y = "SHM Ratio", fill = "Isotype") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(lineheight = .8, size = 20, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 18, colour = "black"),
      axis.text = element_text(face = "bold", size = 14, colour = "black"),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 14)
    ) +
    scale_fill_manual(values = paste("#", c("ED7D31", "FFC000", "954F72", "0563C1"), sep = "")) +
    coord_cartesian(clip = "off") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T,
      step.increase = 0.1
    )
}

SHM.clu <- function(bcr_clusters) {
  ratio.lis <- c()
  for (kk in bcr_clusters) {
    seqs0 <- kk[["junction"]]
    MSAalign <- msa(Biostrings::DNAStringSet(seqs0), "ClustalW")
    seqs <- as.character(attributes(MSAalign)$unmasked)
    msa.l <- unique(nchar(seqs))
    nn <- length(seqs)
    dist <- matrix(0, nn, nn)
    clu.l <- msa.l * nn
    for (ii in 1:nn) {
      for (jj in ii:nn) {
        if (jj == ii) {
          next
        }
        tmp.dist <- SeqDist(seqs[ii], seqs[jj])
        dist[jj, ii] <- dist[jj, ii] + tmp.dist
      }
    }
    dist[dist > 1] <- 0
    dist.is <- sum(dist)
    tmp.ratio <- dist.is / clu.l
    ratio.lis <- c(ratio.lis, tmp.ratio)
  }
  return(ratio.lis)
}

#' Function: Scatter diagram of cluster size (x axis), cluster number (y axis) and SHM ratio (point size) of all clusters in a sample
#'
#' Plot scatter diagram of all clusters in a sample. Each point represents clusters with a given clusters size.
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return A scatter plot of cluster size and cluster number of all clusters inferred from a repertoire. The size of point is positively correlated with the SHM ratio.
#' @export
#'
#' @examples
#' data("bcr_clusters")
#' SHM.cluster(bcr_clusters)
SHM.cluster <- function(bcr_clusters) {
  clu.ratio <- SHM.clu(bcr_clusters)
  clu.size <- sapply(bcr_clusters, function(x) nrow(x))
  clu_size <- as.numeric(names(table(clu.size)))
  clu_n <- as.numeric(table(clu.size))
  clu_ratio <- c()
  for (si in clu_size) {
    tmp.loc <- which(clu.size %in% si)
    tmp.ratio <- mean(clu.ratio[tmp.loc])
    clu_ratio <- c(clu_ratio, tmp.ratio)
  }
  df <- data.frame(
    Cluster_size = clu_size,
    N_Clusters = clu_n,
    SHM_Ratio = clu_ratio
  )
  ggplot(df, aes(x = Cluster_size, y = N_Clusters)) +
    geom_point(aes(size = SHM_Ratio), color = "blue", alpha = 0.5) +
    ylab("Cluster Number") +
    xlab("Cluster Size") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 18, colour = "black"),
      axis.text = element_text(face = "bold", size = 16, colour = "black"),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 16)
    )
}

# CSR sample
#' Function: Network of CSR within a sample
#'
#' Plot CSR network within a sample.
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return Visualization of isotype co-occurrence within BCR clusters. Circle size represents the number of clusters carrying a given isotype. Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses. The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among clusters.
#' @export
#'
#' @examples
#' data("bcr_clusters")
#' CSR.sample(bcr_clusters)
CSR.sample <- function(bcr_clusters) {
  bcr.cluster.isotypes <- NULL
  for (i in seq_along(bcr_clusters)) {
    tmp <- bcr_clusters[[i]]
    tmp.cw <- matrix(0, 1, ncol = 11)
    colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
    tmp.cw[, "Cluster"] <- i
    tmp.is <- tmp[["c_call"]]
    if (length(tmp.is) == 0) tmp.is <- tmp[["C_CALL"]]
    tmp.is <- Iso.fst(tmp.is)
    id <- which(tmp.is %in% c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2"))
    if (length(id) > 0) {
      tmp.is[-id] <- "Unidentified"
    } else {
      tmp.is <- rep("Unidentified", length(tmp.is))
    }
    isotype.count <- table(tmp.is)
    tmp.cw[, names(isotype.count)] <- isotype.count
    bcr.cluster.isotypes <- rbind(bcr.cluster.isotypes, tmp.cw)
  }
  bcr_clusters.count <- length(bcr_clusters)
  network.matrix <- matrix(0, 9, 9)
  rownames(network.matrix) <- colnames(network.matrix) <- c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
  for (i in seq_len(nrow(bcr.cluster.isotypes))) {
    id <- colnames(bcr.cluster.isotypes)[which(as.numeric(bcr.cluster.isotypes[i, 2:10]) > 0) + 1]
    network.matrix[id, id] <- network.matrix[id, id] + 1
  }
  for (i in 1:9) {
    network.matrix[i, i] <- 0
  }
  network.edge <- reshape2::melt(network.matrix)
  network.edge <- network.edge[-which(network.edge$value == 0), ]

  network.data <- network(network.matrix, directed = T, matrix.type = "adjacency")
  network.data %v% "Isotype" <- c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
  set.edge.attribute(network.data, attrname = "as", value = network.edge[, "value"] * 15 / bcr_clusters.count)
  node.attr <- round(rowSums(network.matrix) / sum(network.matrix), 4)
  node.attr[which(node.attr == 0)] <- 0.0001
  node.value <- as.numeric(node.attr)
  cbbPalette2 <- paste("#", c("FFC000", "ED7D31", "0563C1", "44546A", "954F72", "5B9BD5", "E7E6E6", "70AD47", "4472C4"), sep = "")
  GGally::ggnet2(network.data,
    mode = "circle", color = "Isotype", size = "Isotype",
    palette = c("IGHM" = cbbPalette2[1], "IGHD" = cbbPalette2[2], "IGHG3" = cbbPalette2[3], "IGHG1" = cbbPalette2[4], "IGHA1" = cbbPalette2[5], "IGHG2" = cbbPalette2[6], "IGHG4" = cbbPalette2[7], "IGHE" = cbbPalette2[8], "IGHA2" = cbbPalette2[9]),
    size.palette = c("IGHM" = node.value[1], "IGHD" = node.value[2], "IGHG3" = node.value[3], "IGHG1" = node.value[4], "IGHA1" = node.value[5], "IGHG2" = node.value[6], "IGHG4" = node.value[7], "IGHE" = node.value[8], "IGHA2" = node.value[9]),
    edge.size = "as", edge.color = "grey80"
  ) +
    labs(x = "", y = "", color = "Isotype") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      plot.title = element_text(lineheight = .8, size = 20, face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 16)
    ) +
    guides(size = "none")
}

#' Function: Network of CSR within a cluster
#'
#' Plot CSR network within a cluster.
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param n select a cluster to visualize
#'
#' @return Visualization of isotype co-occurrence within a BCR cluster. Circle size represents the number of sequences carrying a given isotype. Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses. The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
#' @export
#'
#' @examples
#' library(fastBCR)
#' data("bcr_clusters")
#' CSR.cluster(bcr_clusters, 1)
CSR.cluster <- function(bcr_clusters, n) {
  tmp <- bcr_clusters[[n]]
  bcr_clusters.count <- nrow(tmp)
  tmp.cw <- matrix(0, 1, ncol = 11)
  colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
  tmp.cw[, "Cluster"] <- n
  tmp.is <- tmp[["c_call"]]
  if (length(tmp.is) == 0) tmp.is <- tmp[["C_CALL"]]
  tmp.is <- Iso.fst(tmp.is)
  id <- which(tmp.is %in% c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2"))
  if (length(id) > 0) {
    tmp.is[-id] <- "Unidentified"
  } else {
    tmp.is <- rep("Unidentified", length(tmp.is))
  }
  isotype.count <- table(tmp.is)
  tmp.cw[, names(isotype.count)] <- isotype.count
  network.matrix <- matrix(0, 9, 9)
  rownames(network.matrix) <- colnames(network.matrix) <- c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
  for (i in seq_len(nrow(tmp.cw))) {
    id <- colnames(tmp.cw)[which(as.numeric(tmp.cw[i, 2:10]) > 0) + 1]
    network.matrix[id, id] <- network.matrix[id, id] + 1
  }
  for (i in 1:9) {
    network.matrix[i, i] <- 0
  }
  network.edge <- reshape2::melt(network.matrix)
  network.edge <- network.edge[-which(network.edge$value == 0), ]
  network.data <- network(network.matrix, directed = T, matrix.type = "adjacency")
  network.data %v% "Isotype" <- c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
  set.edge.attribute(network.data, attrname = "as", value = network.edge[, "value"] * 15 / bcr_clusters.count)
  node.attr <- round(rowSums(network.matrix) / sum(network.matrix), 4)
  node.attr[which(node.attr == 0)] <- 0.0001
  node.value <- as.numeric(node.attr)
  cbbPalette2 <- paste("#", c("FFC000", "ED7D31", "0563C1", "44546A", "954F72", "5B9BD5", "E7E6E6", "70AD47", "4472C4"), sep = "")
  GGally::ggnet2(network.data,
    mode = "circle", color = "Isotype",
    palette = c("IGHM" = cbbPalette2[1], "IGHD" = cbbPalette2[2], "IGHG3" = cbbPalette2[3], "IGHG1" = cbbPalette2[4], "IGHA1" = cbbPalette2[5], "IGHG2" = cbbPalette2[6], "IGHG4" = cbbPalette2[7], "IGHE" = cbbPalette2[8], "IGHA2" = cbbPalette2[9]),
    size = "Isotype", size.palette = c("IGHM" = node.value[1], "IGHD" = node.value[2], "IGHG3" = node.value[3], "IGHG1" = node.value[4], "IGHA1" = node.value[5], "IGHG2" = node.value[6], "IGHG4" = node.value[7], "IGHE" = node.value[8], "IGHA2" = node.value[9]),
    edge.size = "as", edge.color = "grey80") +
    labs(x = "", y = "", color = "Isotype") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(plot.title = element_text(lineheight = .8, size = 20, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 18),
          legend.text = element_text(face = "bold", size = 16)) +
    guides(size = "none")
}
