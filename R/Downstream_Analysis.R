#' Plot: Tree map of clusters
#'
#' @param input BCR repertoire
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return tree map
#' @export
treemap.plot <- function(input, bcr_clusters) {
  size <- sapply(bcr_clusters, function(x) nrow(x))
  clustered_ids <- unique(unlist(lapply(bcr_clusters, function(x) x[["sequence_id"]])))
  unclu_num <- nrow(input) - length(clustered_ids)
  # size <- c(size, rep(1, unclu_num))
  size <- c(size, unclu_num)
  data <- data.frame(
    category = 1:length(size),
    size = size
  )

  treemap::treemap(data,
    index = "category", vSize = "size", vColor = "size",
    type = "size", fontsize.labels = 0, border.col = "gray"
  )
}

#' Function: Merge clustered sequences from all BCR clusters and remove duplicates
#'
#' @param bcr_cluster clonal families inferred by fastBCR
#'
#' @return clustered sequences from all BCR clusters
#' @export
clu2df <- function(bcr_cluster) {
  ss <- length(bcr_cluster)
  cluster.data <- bcr_cluster[[1]][, c("sequence_id", "v_call", "j_call", "junction", "junction_aa")]
  if (ss > 1) {
    for (i in 2:ss) {
      tmp <- bcr_cluster[[i]][, c("sequence_id", "v_call", "j_call", "junction", "junction_aa")]
      cluster.data <- rbind(cluster.data, tmp)
    }
  }
  cluster.data <- cluster.data[!duplicated(cluster.data$sequence_id), ]
  return(cluster.data)
}

#' Calculation: Data.frame for metrics of clusters
#'
#' @param input BCR repertoire
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param group_label Label used to divide different groups
#'
#' @return Data.frame for metrics of bcr clusters, including Tcf20 score, average SHM ratio of clusters and SHM ratios of four isotypes (IGHD, IGHM, IGHA, IGHG).
#' @export
Metrics.df <- function(input, bcr_clusters, group_label) {
  seq_num <- nrow(input)
  clu_num <- length(bcr_clusters)

  # Tcf20 score
  top20 <- bcr_clusters[1:20]
  Tcf20 <- nrow(clu2df(top20)) / seq_num

  # cluster SHM
  SHM_clu <- SHM.clu(bcr_clusters)
  SHM_mean <- mean(SHM_clu)

  # isotype SHM
  SHM_iso <- SHM.iso(bcr_clusters)

  df <- list(
    Group = group_label,
    Tcf20 = round(Tcf20, 3),
    SHM_ratio = round(SHM_mean, 5),
    SHM_iso = SHM_iso
  )

  return(df)
}

#' Plot: Bubble plot between the size and number of clusters between the two groups
#'
#' @param df Data.frame for the size and number of clusters between the two groups
#'
#' @return Bubble plot
#' @export
clu.size.plot <- function(df) {
  ggplot(df, aes(x = size, y = num, fill = group)) +
    geom_point(aes(size = size), alpha = 0.75, shape = 21, color = "black") +
    labs(x = "Size of Clusters", y = "Number of Clusters", fill = "group", size = "Size of Clusters") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_fill_manual(values = c("#ED7D31", "#0563C1")) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(clip = "off")
}

#' Plot: Boxplot of the Tcf20 scores between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param df Data.frame for the Tcf20 score of clusters between the two groups
#'
#' @return Boxplot
#' @export
Tcf20.plot <- function(df) {
  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = df, aes(x = group, y = value, fill = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6), color = "black") +
    geom_boxplot(position = position_dodge(0.6), linewidth = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = group), position = position_jitter(0.1), shape = 21, size = 2, alpha = 0.9) +
    xlab("") +
    ylab("Tcf20 score") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = c("#ED7D31", "#0563C1")) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T
    ) +
    coord_cartesian(clip = "off")
}

#' Plot: Pieplot of V/J gene usage frequency
#'
#' @param clustered_seqs Clustered sequences from a sample
#' @param colname "v_call" for V gene or "j_call" for J gene
#'
#' @return Pie chart. The top ten most frequent genes are shown, and the rest are represented by 'others'
#' @export
pie.freq.plot <- function(clustered_seqs, colname = c("v_call", "j_call")) {
  gene_table <- as.data.frame(table(clustered_seqs[, colname]))
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

  ggplot(tmp_table, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(x = "", y = "", title = colname, fill = "") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(lineheight = .8, size = 20, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_fill_manual(breaks = tmp_table$Var1, values = RColorBrewer::brewer.pal(length(pie.number), "Paired")) +
    guides(fill = guide_legend(byrow = TRUE))
}

#' Plot: Heatmap of VJpair gene usage frequency within a sample
#'
#' @param clustered_seqs Clustered sequences from a sample
#'
#' @return Heatmap
#' @export
vjpair.sample.plot <- function(clustered_seqs) {
  vjpair <- paste(clustered_seqs$v_call, clustered_seqs$j_call, sep = "_")
  clustered_seqs$vjpair <- vjpair
  uni_V <- names(table(clustered_seqs$v_call))
  uni_J <- names(table(clustered_seqs$j_call))
  all_VJ <- c()
  for (i in 1:length(uni_J)) {
    all_VJ <- c(all_VJ, paste(uni_V, uni_J[i], sep = "_"))
  }

  df <- sort(table(vjpair), decreasing = T) %>%
    prop.table() %>%
    as.data.frame() %>%
    rbind(data.frame(
      vjpair = setdiff(all_VJ, vjpair),
      Freq = 0
    )) %>%
    mutate(
      v_call = unlist(strsplit(as.character(vjpair), "_"))[seq(1, 2 * length(vjpair), 2)],
      j_call = unlist(strsplit(as.character(vjpair), "_"))[seq(2, 2 * length(vjpair), 2)]
    ) %>%
    arrange(v_call, j_call) %>%
    select(v_call, j_call, Freq) %>%
    reshape2::dcast(v_call ~ j_call)
  rownames(df) <- df$v_call
  df <- df[, -1]

  colormap <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  breaks <- seq(min(unlist(c(df))), max(unlist(c(df))), length.out = 100)
  pheatmap::pheatmap(df,
    color = colormap,
    breaks = breaks,
    border_color = "black"
  )
}

#' Plot: Heatmap of VJpair gene usage frequency between groups
#'
#' @param group1_clustered_seqs Clustered sequences from group1
#' @param group1_label Label of group1
#' @param group2_clustered_seqs Clustered sequences from group2
#' @param group2_label Label of group2
#'
#' @return Heatmap
#' @export
vjpair.group.plot <- function(group1_clustered_seqs, group1_label, group2_clustered_seqs, group2_label) {
  group1_clustered_seqs$vjpair <- paste(group1_clustered_seqs$v_call, group1_clustered_seqs$j_call, sep = "_")
  group2_clustered_seqs$vjpair <- paste(group2_clustered_seqs$v_call, group2_clustered_seqs$j_call, sep = "_")
  group1_V <- table(group1_clustered_seqs$v_call)
  group2_V <- table(group2_clustered_seqs$v_call)
  group1_J <- table(group1_clustered_seqs$j_call)
  group2_J <- table(group2_clustered_seqs$j_call)
  uni_V <- union(names(group1_V), names(group2_V))
  uni_J <- union(names(group1_J), names(group2_J))
  all_VJ <- c()
  for (i in 1:length(uni_J)) {
    all_VJ <- c(all_VJ, paste(uni_V, uni_J[i], sep = "_"))
  }

  df1 <- sort(table(group1_clustered_seqs$vjpair), decreasing = T) %>%
    prop.table() %>%
    as.data.frame() %>%
    rbind(data.frame(
      Var1 = setdiff(all_VJ, group1_clustered_seqs$vjpair),
      Freq = 0
    )) %>%
    mutate(
      v_call = unlist(strsplit(as.character(Var1), "_"))[seq(1, 2 * length(all_VJ), 2)],
      j_call = unlist(strsplit(as.character(Var1), "_"))[seq(2, 2 * length(all_VJ), 2)]
    ) %>%
    arrange(v_call, j_call) %>%
    select(v_call, j_call, Freq) %>%
    reshape2::dcast(v_call ~ j_call)
  rownames(df1) <- df1$v_call
  df1 <- df1[, -1]

  df2 <- sort(table(group2_clustered_seqs$vjpair), decreasing = T) %>%
    prop.table() %>%
    as.data.frame() %>%
    rbind(data.frame(
      Var1 = setdiff(all_VJ, group2_clustered_seqs$vjpair),
      Freq = 0
    )) %>%
    mutate(
      v_call = unlist(strsplit(as.character(Var1), "_"))[seq(1, 2 * length(all_VJ), 2)],
      j_call = unlist(strsplit(as.character(Var1), "_"))[seq(2, 2 * length(all_VJ), 2)]
    ) %>%
    arrange(v_call, j_call) %>%
    select(v_call, j_call, Freq) %>%
    reshape2::dcast(v_call ~ j_call)
  rownames(df2) <- df2$v_call
  df2 <- df2[, -1]

  colormap <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  breaks <- seq(min(unlist(c(df1, df2))), max(unlist(c(df1, df2))), length.out = 100)
  p1 <- pheatmap::pheatmap(df1,
    color = colormap,
    breaks = breaks,
    border_color = "black",
    legend = TRUE
  )
  p2 <- pheatmap::pheatmap(df2,
    color = colormap,
    breaks = breaks,
    border_color = "black",
    legend = TRUE
  )
  cowplot::plot_grid(p1$gtable, p2$gtable, align = "vh", labels = c(group1_label, group2_label), ncol = 2)
}

#' Calculation: Gene usage frequency of clustered sequences from a sample
#'
#' @param clustered_seqs Clustered sequences from a sample
#' @param colname "v_call" for V gene or "j_call" for J gene
#' @param uni_gene All unique genes in two groups of samples
#' @param group_label Label used to divide different groups
#'
#' @return Data.frame for gene usage frequency
#' @export
gene.fre.df <- function(clustered_seqs, colname = c("v_call", "j_call"), uni_gene, group_label) {
  df <- table(clustered_seqs[, colname]) %>%
    prop.table() %>%
    as.data.frame()
  if (length(setdiff(uni_gene, clustered_seqs[, colname])) != 0) {
    df <- rbind(df, data.frame(
      Var1 = setdiff(uni_gene, clustered_seqs[, colname]),
      Freq = 0
    ))
  }
  df <- mutate(df, x = as.character(Var1)) %>%
    arrange(x) %>%
    dplyr::rename(y = Freq) %>%
    mutate(group = group_label) %>%
    select(-Var1) %>%
    select(x, everything())

  return(df)
}

#' Plot: Boxplot of the gene usage frequency between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param df Data.frame for the gene usage frequency between the two groups
#'
#' @return Boxplot
#' @export
gene.fre.plot <- function(df) {
  df$y <- df$y * 100
  df$x <- as.factor(df$x)
  stat.test <- df %>%
    group_by(x) %>%
    rstatix::wilcox_test(y ~ group) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x = "x", dodge = 1)
  ggpubr::ggboxplot(df,
    x = "x", y = "y", fill = "group",
    outlier.shape = NA, palette = c("#ED7D31", "#0563C1")
  ) +
    ggpubr::stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0, hide.ns = TRUE) +
    labs(x = "", y = "Frequency (%)", fill = "Group") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    coord_cartesian(clip = "off")
}

#' Plot: Histogram and density plot of junction amino acid length of clustered sequences from a sample
#'
#' @param clustered_seqs Clustered sequences from a sample
#'
#' @return junction amino acid length distribution from a sample
#' @export
len.sample.plot <- function(clustered_seqs) {
  nchar_junc <- nchar(clustered_seqs[, "junction_aa"])
  junc <- as.clustered_seqs.frame(nchar_junc)

  ggplot(junc, aes(x = nchar_junc)) +
    geom_histogram(
      binwidth = 2, colour = "black",
      aes(y = after_stat(density), fill = after_stat(count))
    ) +
    stat_function(fun = dnorm, color = "darkblue", size = 2, args = list(mean = mean(nchar_junc), sd = sd(nchar_junc))) +
    labs(x = "Length", y = "Density") +
    ggtitle(paste("Junction Length Distribution", sep = "")) +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    )
}

#' Plot: Density ridges of junction amino acid length
#'
#' @param df Data.frame for the junction amino acid length of clustered sequences between the two groups
#'
#' @return junction amino acid length distribution between the two groups
#' @export
len.group.plot <- function(df) {
  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(df, aes(x = value, y = group, fill = group)) +
    ggridges::geom_density_ridges(aes(fill = group),
      colour = "black", alpha = 0.8, quantile_lines = TRUE,
      vline_linetype = "dashed", quantiles = 2, bandwidth = 1.5
    ) +
    labs(x = "length", y = "") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = c("#ED7D31", "#0563C1")) +
    scale_x_continuous(
      limits = c(5, 30),
      breaks = c(5, 10, 15, 20, 25, 30)
    ) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      bracket.size = 1,
      hide.ns = T
    ) +
    coord_cartesian(clip = "off")
}

name <- function(bcr_clusters, index, loc) {
  id <- bcr_clusters[[index]][["sequence_id"]]
  id <- id[loc]
  c_call <- bcr_clusters[[index]][["c_call"]]
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

#' Plot: Visualization of multiple sequence alignment (MSA) of junction sequences within a cluster
#'
#' Perform MSA on the junction of all sequences in the clonal family for visualization
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param index index of cluster
#' @param type 'DNA' or 'AA'
#'
#' @return Visualization of MSA. From top to bottom are title("Vgene_Jgene_Length"), seqlogo, ggmsa, and msabar
#' @export
msa.plot <- function(bcr_clusters, index, type = c("AA", "DNA")) {
  if (type == "AA") {
    seqs_aa0 <- bcr_clusters[[index]][["junction_aa"]]
    MSAalign_aa <- msa::msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW")
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::AAStringSet(seqs_aa)
  } else if (type == "DNA") {
    seqs_aa0 <- bcr_clusters[[index]][["junction"]]
    MSAalign_aa <- msa::msa(Biostrings::DNAStringSet(seqs_aa0), "ClustalW")
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::DNAStringSet(seqs_aa)
  }

  seqs_aa_no <- gsub("-", "", seqs_aa)
  loc <- sort.msa(seqs_aa_no, seqs_aa0)
  SEQs_aa@ranges@NAMES <- name(bcr_clusters, index, loc)
  v_call <- unique(bcr_clusters[[index]][["v_call"]])
  j_call <- unique(bcr_clusters[[index]][["j_call"]])
  len <- unique(nchar(seqs_aa))
  vjlen <- paste(paste(v_call, collapse = ","), "_", paste(j_call, collapse = ","), "_Length:", len, sep = "")
  clu.title <- vjlen
  ggmsa::ggmsa(SEQs_aa, char_width = 0.6, seq_name = TRUE) +
    ggtitle(label = clu.title) +
    theme(
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0),
      axis.title = element_text(size = 14, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(hjust = 0)
    ) +
    geom_seqlogo(adaptive = FALSE) +
    geom_msaBar()
}

#' Plot: Visualization of seqlogo of junction sequences within a cluster
#'
#' Perform seqlogo on the junction of all sequences in the clonal family for visualization
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param index index of cluster
#' @param type 'DNA' or 'AA'
#'
#' @return Seqlogo
#' @export
seqlogo.plot <- function(bcr_clusters, index, type = c("AA", "DNA")) {
  if (type == "AA") {
    seqs_aa0 <- bcr_clusters[[index]][["junction_aa"]]
    MSAalign_aa <- msa::msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW")
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::AAStringSet(seqs_aa)
  } else if (type == "DNA") {
    seqs_aa0 <- bcr_clusters[[index]][["junction"]]
    MSAalign_aa <- msa::msa(Biostrings::DNAStringSet(seqs_aa0), "ClustalW")
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::DNAStringSet(seqs_aa)
  }

  seqs_aa_no <- gsub("-", "", seqs_aa)
  loc <- sort.msa(seqs_aa_no, seqs_aa0)
  SEQs_aa@ranges@NAMES <- name(bcr_clusters, index, loc)
  v_call <- unique(bcr_clusters[[index]][["v_call"]])
  j_call <- unique(bcr_clusters[[index]][["j_call"]])
  len <- unique(nchar(seqs_aa))
  vjlen <- paste(paste(v_call, collapse = ","), "_", paste(j_call, collapse = ","), "_Length:", len, sep = "")
  clu.title <- vjlen
  ggseqlogo::ggseqlogo(seqs_aa, method = "prob") +
    ggtitle(clu.title) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0.5),
      strip.text = element_text(lineheight = .8, size = 12, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    coord_cartesian(clip = "off")
}

#' Plot: Reconstructing B cell lineage trees with minimum spanning tree and genotype abundances (ClonalTree)
#'
#' @param  bcr_clusters clonal families inferred by fastBCR
#' @param index index of cluster
#'
#' @return the reconstructed BCR lineage tree
#' @export
clonal.tree.plot <- function(bcr_clusters, index) {
  ids <- bcr_clusters[[index]]$sequence_id
  l <- nchar(ids)
  if (min(l) > 5) {
    ids <- stringr::str_sub(ids, l - 4, l)
  }
  seqs <- bcr_clusters[[index]]$junction
  MSAalign_dna <- msa::msa(Biostrings::DNAStringSet(seqs), "ClustalW")
  seqs_dna <- as.character(attributes(MSAalign_dna)$unmasked)
  splt <- strsplit(seqs_dna, "")
  DNA_len <- length(splt[[1]])
  most_DNA <- c()
  for (j in 1:DNA_len) {
    DNA_lis <- sapply(splt, function(x) x[j])
    tmp.DNA <- names(sort(table(DNA_lis), decreasing = TRUE))[1]
    if (tmp.DNA == "-") {
      tmp.DNA <- names(sort(table(DNA_lis), decreasing = TRUE))[2]
    }
    most_DNA <- c(most_DNA, tmp.DNA)
  }
  consensus_DNA <- paste(most_DNA, collapse = "")
  con.loc <- which(consensus_DNA %in% seqs)

  if (length(con.loc) == 0) {
    ids <- c("naive", ids)
    seqs <- c(consensus_DNA, seqs)
  } else {
    ids[con.loc] <- "naive"
  }

  fasta <- paste0(">", ids, "\n", seqs)
  write.table(fasta,
    file = "ClonalTree/Examples/input/CF.fasta",
    row.names = F,
    col.names = F,
    quote = F
  )
  system("python3 ClonalTree/src/clonalTree.py  -i ClonalTree/Examples/input/CF.fasta -o ClonalTree/Examples/output/CF.abRT.nk -a 1 -r 1 -t 1")
  tree <- ape::read.tree("ClonalTree/Examples/output/CF.abRT.nk")
  ggtree::ggtree(tree, size = 1, col = "deepskyblue3", options(ignore.negative.edge = TRUE)) +
    ggtree::geom_nodelab(size = 4, color = "purple4", hjust = 1.2, vjust = -0.3) +
    ggtree::geom_nodepoint(size = 3, color = "orange2") +
    ggtree::geom_tiplab(size = 4, color = "purple3", hjust = -0.3) +
    ggtree::geom_tippoint(size = 2, color = "deepskyblue3") +
    ggtree::theme_tree2() +
    theme(axis.text = element_text(face = "bold", size = 12, colour = "black"))
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

#' Calculation: SHM ratios of clusters from a sample
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return SHM ratios of clusters
#' @export
SHM.clu <- function(bcr_clusters) {
  ratio.lis <- c()
  for (kk in bcr_clusters) {
    seqs0 <- kk[["junction"]]
    MSAalign <- msa::msa(Biostrings::DNAStringSet(seqs0), "ClustalW")
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

#' Plot: Boxplot of the SHM ratios between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param df Data.frame for the SHM ratios of clusters between the two groups
#'
#' @return Boxplot
#' @export
SHM.plot <- function(df) {
  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = df, aes(x = group, y = value, fill = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6), color = "black") +
    geom_boxplot(position = position_dodge(0.6), linewidth = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = group), position = position_jitter(0.1), shape = 21, size = 2, alpha = 0.9) +
    xlab("") +
    ylab("SHM ratio") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = c("#ED7D31", "#0563C1")) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T
    ) +
    coord_cartesian(clip = "off")
}

#' Calculation: SHM ratios of clustered sequences in different isotypes from a sample
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return SHM ratios of clustered sequences in four isotypes ("IGHD", "IGHM", "IGHA", and "IGHG")
#' @export
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
      MSAalign <- msa::msa(Biostrings::DNAStringSet(tmp.seqs0), "ClustalW")
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

#' Plot: Boxplot of SHM ratios of clustered sequences in different isotypes from a group
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param df Data.frame for SHM ratios of clustered sequences in four isotypes ("IGHD", "IGHM", "IGHA", and "IGHG") from multiple samples
#'
#' @return Boxplot
#' @export
SHM.iso.plot <- function(df) {
  res <- ggpubr::compare_means(formula = value ~ Isotypes, df) %>% filter(p.signif != "ns")
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
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = paste("#", c("ED7D31", "FFC000", "954F72", "0563C1"), sep = "")) +
    coord_cartesian(clip = "off") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = TRUE,
      step.increase = 0.1
    )
}

#' Plot: Network of CSR within a cluster
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param index index of cluster
#'
#' @return Visualization of isotype co-occurrence within a BCR cluster. Circle size represents the number of sequences carrying a given isotype. Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses. The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
#' @export
CSR.cluster.plot <- function(bcr_clusters, index) {
  tmp <- bcr_clusters[[index]]
  bcr_clusters.count <- nrow(tmp)
  tmp.cw <- matrix(0, 1, ncol = 11)
  colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
  tmp.cw[, "Cluster"] <- index
  tmp.is <- tmp[["c_call"]]
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
  network.data <- network(network.matrix, directed = TRUE, matrix.type = "adjacency")
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
    edge.size = "as", edge.color = "grey80"
  ) +
    labs(x = "", y = "", color = "Isotype") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    guides(size = "none")
}

#' Plot: Network of CSR within a sample
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#'
#' @return Visualization of isotype co-occurrence within BCR clusters. Circle size represents the number of clusters carrying a given isotype. Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses. The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among clusters.
#' @export
CSR.sample.plot <- function(bcr_clusters) {
  bcr.cluster.isotypes <- NULL
  for (i in seq_along(bcr_clusters)) {
    tmp <- bcr_clusters[[i]]
    tmp.cw <- matrix(0, 1, ncol = 11)
    colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
    tmp.cw[, "Cluster"] <- i
    tmp.is <- tmp[["c_call"]]
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

  network.data <- network(network.matrix, directed = TRUE, matrix.type = "adjacency")
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
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    guides(size = "none")
}

findNa <- function(sampledata, NAb_v, NAb_j, NAb_cdr3) {
  nn <- nrow(sampledata)
  ids <- c()
  for (i in 1:nn) {
    sample_id <- sampledata$sequence_id[i]
    sample_v <- sampledata$v_call[i]
    sample_v <- unlist(strsplit(sample_v, "\\*"))[1]
    sample_j <- sampledata$j_call[i]
    sample_j <- unlist(strsplit(sample_j, "\\*"))[1]
    NAb_loc <- which(NAb_v == sample_v & NAb_j == sample_j)
    if (length(NAb_loc) == 0) {
      next
    }
    tmp_NAb_cdr3 <- NAb_cdr3[NAb_loc]
    sample_cdr3 <- sampledata$junction_aa[i]
    cdr3_l <- nchar(sample_cdr3)
    sample_cdr3 <- substr(sample_cdr3, 2, cdr3_l - 1)
    if (length(grep(paste("^", sample_cdr3, "$", sep = ""), tmp_NAb_cdr3)) != 0) {
      ids <- c(ids, sample_id)
    }
  }
  ids <- unique(ids)

  return(ids)
}

#' Calculation: NAb ratio (the ratio derived from dividing the tally of NAb sequences encompassed by clonal families by the aggregate count of NAb sequences present within each sample)
#'
#' @param input BCR repertoire
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param NAb_v V gene of all NAb sequences
#' @param NAb_j J gene of all NAb sequences
#' @param NAb_cdr3 CDR3 sequences of all NAb sequences
#'
#' @return NAb ratio
#' @export
NAb.ratio <- function(input, bcr_clusters, NAb_v, NAb_j, NAb_cdr3) {
  clu_data <- clu2df(bcr_clusters)
  id_seq <- findNa(input, NAb_v, NAb_j, NAb_cdr3)
  id_clu <- findNa(clu_data, NAb_v, NAb_j, NAb_cdr3)
  if (length(id_seq) == 0) {
    NAb_ratio <- 0
  } else {
    NAb_ratio <- length(id_clu) / length(id_seq)
  }

  return(NAb_ratio)
}

#' Plot: Boxplot of the NAb ratios between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param df Data.frame for NAb ratios of samples
#'
#' @return Boxplot
#' @export
NAb.ratio.plot <- function(df) {
  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = df, aes(x = group, y = value, fill = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5, position = position_dodge(0.6), color = "black") +
    geom_boxplot(position = position_dodge(0.6), size = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = group), position = position_jitter(0.1), shape = 21, size = 2, alpha = 0.9) +
    xlab("") +
    ylab("NAb Ratio") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = c("#ED7D31", "#0563C1")) +
    coord_cartesian(clip = "off") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T,
      step.increase = 0.1
    ) +
    scale_y_continuous(limits = c(0, 1))
}

#' Plot: ROC curve of the NAb ratios between the two groups.
#'
#' @param roc a “roc” object, a list of class “roc” buliding from pROC package
#'
#' @return ROC curve
#' @export
NAb.roc.plot <- function(roc) {
  lab <- paste("AUC: ", round(roc$auc, 3), sep = "")
  pROC::ggroc(roc, size = 1, legacy.axes = T) +
    theme_classic() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
      colour = "grey",
      linetype = "dotdash"
    ) +
    annotate("text", x = 0.2, y = 1, label = lab, size = 5, color = "black") +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    ) +
    coord_cartesian(clip = "off")
}
