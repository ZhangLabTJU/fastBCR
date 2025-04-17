#' Function: Classification of clustered and unclustered sequences
#'
#' @param pro_data_list A list where each element is the preprocessed data named after its filename
#' @param clusters_list A list where each element is a list of clonal families inferred by fastBCR from a single sample
#'
#' @return A list containing clustered and unclustered sequences in each sample
#' 'clustered_seqs' is a list where each element is the clustered sequences in a sample named after its filename
#' 'unclustered_seqs' is a list where each element is the unclustered sequences in a sample named after its filename
#' @export
Clustered.seqs <- function(pro_data_list, clusters_list) {
  clustered_seqs_sample <- list()
  unclustered_seqs_sample <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(clusters_list), style = 3)
  for (i in seq_along(clusters_list)) {
    var_name <- names(clusters_list[i])
    pro_data <- pro_data_list[[i]]
    clusters <- clusters_list[[i]]
    clustered_index <- unique(unlist(sapply(clusters, function(x) x$clonotype_index)))
    clustered_seqs <- pro_data[clustered_index, ]
    all_index <- pro_data$clonotype_index
    unclustered_seqs <- pro_data[-clustered_index, ]
    clustered_seqs_sample[[var_name]] <- clustered_seqs
    unclustered_seqs_sample[[var_name]] <- unclustered_seqs
    setTxtProgressBar(pb, i)
  }
  clustered_seqs <- list(
    clustered_seqs = clustered_seqs_sample,
    unclustered_seqs = unclustered_seqs_sample
  )

  cat(paste0(
    paste0("\nYou have classified clustered and unclustered sequences for ", length(pro_data_list), " samples\n"),
    paste(paste0(
      "Sample '", names(pro_data_list), "' contains ",
      sapply(clustered_seqs$clustered_seqs, function(x) nrow(x)), " clustered sequences and ",
      sapply(clustered_seqs$unclustered_seqs, function(x) nrow(x)), " unclustered sequences"
    ), collapse = "\n")
  ))

  return(clustered_seqs)
}

#' Function: Summarize the number of clusters, the average size of clusters and the proportion of clustered sequences
#'
#' @param pro_data_list A list where each element is the preprocessed data named after its filename
#' @param clusters_list A list where each element is a list of clonal families inferred by fastBCR from a single sample
#'
#' @return A list where each element is the summary of clusters from a single sample
#' @export
Clusters.summary <- function(pro_data_list, clusters_list) {
  Clusters_summary_list <- list()
  for (i in seq_along(pro_data_list)) {
    pro_data <- pro_data_list[[i]]
    pro_data_n <- nrow(pro_data)
    cluster_data <- clusters_list[[i]]
    clustered_data <- unique(unlist(sapply(cluster_data, function(x) x$clonotype_index)))
    clustered_data_n <- length(clustered_data)
    var_name <- names(pro_data_list[i])
    cluster_size <- sapply(cluster_data, function(x) nrow(x))
    Occupying <- scales::percent(clustered_data_n / pro_data_n, accuracy = 0.01)
    clusters_num <- length(cluster_size)
    average_cluster_size <- sum(cluster_size) / clusters_num
    Clusters_summary <- NULL
    Clusters_summary[["number of clusters"]] <- clusters_num
    Clusters_summary[["average size of clusters"]] <- round(average_cluster_size, 2)
    Clusters_summary[["number of clustered seqs"]] <- clustered_data_n
    Clusters_summary[["number of all seqs"]] <- pro_data_n
    Clusters_summary[["proportion of clustered sequences"]] <- Occupying
    Clusters_summary_list[[var_name]] <- Clusters_summary
  }

  return(Clusters_summary_list)
}

#' Function: Calculate the size and number of clusters for the samples within each group
#'
#' @param clusters_list A list where each element of the list is the clonal families inferred by fastBCR for each sample
#' @param group_label Label of group
#'
#' @return Dataframe containing information on the group, cluster size, and number
#' @export
clu.size.df <- function(clusters_list, group_label) {
  clusters_all <- numeric(0)
  for (i in seq_along(clusters_list)) {
    cluster_data <- clusters_list[[i]]
    clusters_all <- c(clusters_all, cluster_data)
  }
  clu.size <- sapply(clusters_all, function(x) nrow(x))
  clu_size <- as.numeric(names(table(clu.size)))
  clu_n <- as.numeric(table(clu.size))
  clu_size_df <- data.frame(
    group = rep(group_label, length(clu_size)),
    size = clu_size,
    num = clu_n
  )

  return(clu_size_df)
}

#' Plot: Bubble plot between the size and number of clusters between the two groups
#'
#' @param clusters_list1 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group1
#' @param group1_label Label of group1
#' @param clusters_list2 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group2
#' @param group2_label Label of group2
#'
#' @return Bubble plot showing the size and number of clusters between the two groups
#' @export
clu.size.plot <- function(clusters_list1, group1_label,
                          clusters_list2, group2_label) {
  df1 <- clu.size.df(clusters_list1, group1_label)
  df2 <- clu.size.df(clusters_list2, group2_label)
  df <- rbind(df1, df2)

  res <- ggpubr::compare_means(formula = size ~ group, df)
  print(paste0("p-value: ", res$p.format))

  ggplot(df, aes(x = size, y = num, fill = group, color = group)) +
    geom_point(aes(size = size), alpha = 0.75, shape = 21) +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(clip = "off")
}

Tcf20.df <- function(clusters_list, group_label) {
  Tcf20_values <- c()
  name <- c()
  for (i in seq_along(clusters_list)) {
    clusters <- clusters_list[[i]]
    top20 <- clusters[1:20]
    top20_n <- sum(sapply(top20, function(x) nrow(x)))
    all_n <- sum(sapply(clusters, function(x) nrow(x)))
    Tcf20 <- top20_n / all_n
    var_name <- names(clusters_list[i])
    name <- c(name, var_name)
    Tcf20_values <- c(Tcf20_values, Tcf20)
  }
  Tcf20_df <- data.frame(
    group = rep(group_label, length(clusters_list)),
    Tcf20 = Tcf20_values
  )
  rownames(Tcf20_df) <- name

  return(Tcf20_df)
}

#' Plot: Boxplot of the Tcf20 scores between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param clusters_list1 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group1
#' @param group1_label Label of group1
#' @param clusters_list2 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group2
#' @param group2_label Label of group2
#'
#' @return Boxplot showing the Tcf20 scores between the two groups
#' @export
Tcf20.plot <- function(clusters_list1, group1_label,
                       clusters_list2, group2_label) {
  df1 <- Tcf20.df(clusters_list1, group1_label)
  df2 <- Tcf20.df(clusters_list2, group2_label)
  Tcf20_df <- rbind(df1, df2)
  print(Tcf20_df)

  res <- ggpubr::compare_means(formula = Tcf20 ~ group, Tcf20_df) %>% filter(p.signif != "ns")
  print(paste0("p-value: ", res$p.format))
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = Tcf20_df, aes(x = group, y = Tcf20, fill = group, color = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6)) +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T
    ) +
    coord_cartesian(clip = "off")
}

#'  Plot: Visualization of clusters from a sample
#'
#' @param pro_data_list A list where each element is the preprocessed data named after its filename
#' @param clusters_list A list where each element is a list of clonal families inferred by fastBCR from a single sample
#' @param index The Index of sample for visualization
#'
#' @return  Point diagram of clusters from a sample where a circle represents a cluster. The size and color of the circle represents the size of the cluster.
#' @export
Clusters.visualization <- function(pro_data_list, clusters_list, index) {
  pro_data <- pro_data_list[[index]]
  cluster_data <- clusters_list[[index]]
  pro_data_name <- names(pro_data_list[index])
  clusters_name <- names(clusters_list[index])
  cluster_size <- sapply(cluster_data, function(x) nrow(x))
  Occupying <- scales::percent(sum(cluster_size) / nrow(pro_data), accuracy = 0.01)
  Occupying <- paste0(Occupying, " (", sum(cluster_size), "/", nrow(pro_data), ")")
  clusters_num <- length(cluster_size)
  set.seed(123)
  position <- ggraph::pack_circles(cluster_size)
  data <- data.frame(
    x = position[, 1],
    y = position[, 2],
    r = sqrt(cluster_size / pi),
    value = cluster_size,
    label = cluster_size
  )
  ggplot(data) +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r, fill = value), data = data, colour = "white") +
    scale_fill_distiller(palette = "Spectral", name = "size") +
    guides(size = FALSE) +
    coord_fixed() +
    theme(panel.border = element_blank()) +
    theme_void() +
    labs(
      title = paste0("Visualization of ", clusters_name, "_clusters (n = ", clusters_num, ")\nOccupying ", Occupying, " of ", pro_data_name)
    ) +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(lineheight = .8, size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(face = "bold", size = 12),
      legend.title = element_text(face = "bold", size = 14)
    )
}

#' Plot: Pieplot of V/J gene usage frequency
#'
#' @param clustered_seqs Clustered sequences
#' @param colname 'v_call' for V gene or 'j_call' for J gene
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
      plot.title = element_text(lineheight = .8, size = 24, face = "bold", hjust = 0.5),
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

#' Plot: Boxplot of the gene usage frequency between the two groups.
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param group1_seqs_list A list containing clustered and unclustered sequences in each sample in group1
#' @param group1_all_clustered_seqs All the clustered sequences in group1
#' @param group1_label Label of group1
#' @param group2_seqs_list A list containing clustered and unclustered sequences in each sample in group2
#' @param group2_all_clustered_seqs All the clustered sequences in group2
#' @param group2_label Label of group2
#' @param colname 'v_call' for V gene or 'j_call' for J gene
#'
#' @return Boxplot showing the V/J gene usage of the clustered sequences between the two groups.
#' @export
gene.fre.plot <- function(group1_seqs_list, group1_all_clustered_seqs, group1_label,
                          group2_seqs_list, group2_all_clustered_seqs, group2_label, colname) {
  group1_seqs_sample <- group1_seqs_list$clustered_seqs
  group2_seqs_sample <- group2_seqs_list$clustered_seqs
  group1_gene <- table(group1_all_clustered_seqs[[colname]])
  group2_gene <- table(group2_all_clustered_seqs[[colname]])
  uni_gene <- union(names(group1_gene), names(group2_gene))
  group1_fre <- Clusters.gene.fre(group1_seqs_sample, colname, uni_gene, group1_label)
  group2_fre <- Clusters.gene.fre(group2_seqs_sample, colname, uni_gene, group2_label)

  gene_fre <- rbind(group1_fre, group2_fre)

  gene_fre$y <- gene_fre$y * 100
  gene_fre$x <- as.factor(gene_fre$x)
  stat.test <- gene_fre %>%
    group_by(x) %>%
    rstatix::wilcox_test(y ~ group) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x = "x", dodge = 1)
  ggpubr::ggboxplot(gene_fre,
    x = "x", y = "y", fill = "group", color = "group",
    outlier.shape = NA
  ) +
    ggpubr::stat_pvalue_manual(stat.test, label = "{p.adj.signif}", tip.length = 0, hide.ns = TRUE) +
    labs(x = "", y = "Frequency (%)") +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    coord_cartesian(clip = "off")
}

#' Plot: Heatmap of V-J gene pair frequency of clustered sequences from a sample
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


#' Plot: Heatmap of the fold change of V-J gene pair frequency between groups
#'
#' Log fold change (log FC) is calculated as the log2 ratio of the average values between group1 and group2 samples
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#' FDR correction was performed with the Benjamini–Hochberg procedure
#'
#' @param group1_seqs_list A list containing clustered and unclustered sequences in each sample in group1
#' @param group1_all_clustered_seqs All the clustered sequences in group1
#' @param group1_label Label of group1
#' @param group2_seqs_list A list containing clustered and unclustered sequences in each sample in group2
#' @param group2_all_clustered_seqs All the clustered sequences in group2
#' @param group2_label Label of group2
#'
#' @return Heatmap showing the fold change of V-J gene pair frequency of clustered sequences between two groups.
#' @export
vjpair.group.plot <- function(group1_seqs_list, group1_all_clustered_seqs, group1_label,
                              group2_seqs_list, group2_all_clustered_seqs, group2_label) {
  group1_seqs_sample <- group1_seqs_list$clustered_seqs
  for (sample in names(group1_seqs_sample)) {
    group1_seqs_sample[[sample]]$vjpair <- paste(group1_seqs_sample[[sample]]$v_call, group1_seqs_sample[[sample]]$j_call, sep = "_")
  }
  group2_seqs_sample <- group2_seqs_list$clustered_seqs
  for (sample in names(group2_seqs_sample)) {
    group2_seqs_sample[[sample]]$vjpair <- paste(group2_seqs_sample[[sample]]$v_call, group2_seqs_sample[[sample]]$j_call, sep = "_")
  }
  group1_all_clustered_seqs$vjpair <- paste(group1_all_clustered_seqs$v_call, group1_all_clustered_seqs$j_call, sep = "_")
  group2_all_clustered_seqs$vjpair <- paste(group2_all_clustered_seqs$v_call, group2_all_clustered_seqs$j_call, sep = "_")
  group1_V <- table(group1_all_clustered_seqs$v_call)
  group2_V <- table(group2_all_clustered_seqs$v_call)
  group1_J <- table(group1_all_clustered_seqs$j_call)
  group2_J <- table(group2_all_clustered_seqs$j_call)
  uni_V <- union(names(group1_V), names(group2_V))
  uni_J <- union(names(group1_J), names(group2_J))
  all_VJ <- c()
  for (i in 1:length(uni_J)) {
    all_VJ <- c(all_VJ, paste(uni_V, uni_J[i], sep = "_"))
  }

  group1_fre <- Clusters.gene.fre(group1_seqs_sample, "vjpair", all_VJ, group1_label)
  group2_fre <- Clusters.gene.fre(group2_seqs_sample, "vjpair", all_VJ, group2_label)

  group1_fre$y[which(group1_fre$y <= 0.001)] <- 0.001
  group2_fre$y[which(group2_fre$y <= 0.001)] <- 0.001

  vjpair.df <- NULL
  for (vj in all_VJ) {
    group1.vj <- group1_fre[which(group1_fre$x == vj), ]
    group1.vj.mean <- mean(group1.vj$y)
    group2.vj <- group2_fre[which(group2_fre$x == vj), ]
    group2.vj.mean <- mean(group2.vj$y)
    vj.df <- rbind(
      group1.vj,
      group2.vj
    )
    res <- ggpubr::compare_means(formula = y ~ group, data = vj.df, p.adjust.method = "BH")
    if (nrow(res) == 0) {
      fdr <- 1
      signif <- "ns"
    } else {
      fdr <- res$p.adj
      signif <- res$p.signif
    }
    tmp.vjpair.df <- data.frame(
      vjpair = vj,
      v_gene = unlist(strsplit(vj, "_"))[1],
      j_gene = unlist(strsplit(vj, "_"))[2],
      group1_mean_value = group1.vj.mean,
      group2_mean_value = group2.vj.mean,
      log = log2(group1.vj.mean / group2.vj.mean),
      FDR = 1 / fdr,
      p.signif = signif
    )
    vjpair.df <- rbind(vjpair.df, tmp.vjpair.df)
  }

  v_gene_n <- length(unique(vjpair.df$v_gene))
  j_gene_n <- length(unique(vjpair.df$j_gene))
  x_seg <- data.frame(
    x = seq(1.5, (v_gene_n - 0.5), 1),
    xend = seq(1.5, (v_gene_n - 0.5), 1),
    y = -Inf,
    yend = Inf
  )
  y_seg <- data.frame(
    x = -Inf,
    xend = Inf,
    y = seq(1.5, (j_gene_n - 0.5), 1),
    yend = seq(1.5, (j_gene_n - 0.5), 1)
  )

  colormap <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  ggplot(vjpair.df, aes(x = v_gene, y = j_gene)) +
    geom_point(aes(size = FDR, color = log), shape = 15) +
    theme_bw() +
    ggtitle(paste(group1_label, group2_label, sep = " / ")) +
    labs(x = NULL, y = NULL, color = expression(log[2] * "F" * "C")) +
    theme(
      panel.grid.major = element_line(colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      plot.title = element_text(lineheight = .8, size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12, colour = "black"),
      axis.text = element_text(face = "bold", size = 10, colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(face = "bold", size = 10)
    ) +
    geom_segment(
      data = x_seg, aes(x = x, xend = xend, y = y, yend = yend),
      color = "black"
    ) +
    geom_segment(
      data = y_seg, aes(x = x, xend = xend, y = y, yend = yend),
      color = "black"
    ) +
    scale_color_gradient2(low = colormap[1], high = colormap[100], mid = "lightgrey") +
    scale_size(breaks = c(1, 20, 100, 1000, 10000), labels = c("ns", "*", "**", "***", "****")) +
    coord_cartesian(clip = "off")
}

#' Calculation: Gene usage frequency of clustered sequences from a sample
#'
#' @param clustered_seqs Clustered sequences from a sample
#' @param colname 'v_call' for V gene or 'j_call' for J gene
#' @param uni_gene All unique genes in two groups of samples
#' @param group_label Label of group
#'
#' @return Data.frame for gene usage frequency
#' @export
gene.fre.df <- function(clustered_seqs, colname, uni_gene, group_label) {
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

#' Function: Calculate the frequency of gene usage in clustered sequences from the samples
#'
#' @param clustered_seqs_sample A list where each element is the clustered sequences for each sample
#' @param colname 'v_call' for V gene or 'j_call' for J gene
#' @param uni_gene 	All unique genes in two groups of samples
#' @param group_label Label of group
#'
#' @return Dataframe containing the gene usage frequency for samples
#' @export
Clusters.gene.fre <- function(clustered_seqs_sample, colname, uni_gene, group_label) {
  clustered_gene_fre <- data.frame()
  for (i in seq_along(clustered_seqs_sample)) {
    clustered_seqs_data <- clustered_seqs_sample[[i]]
    processed_data <- gene.fre.df(clustered_seqs_data, colname, uni_gene, group_label)
    clustered_gene_fre <- rbind(clustered_gene_fre, processed_data)
  }

  return(clustered_gene_fre)
}

#' Plot: Histogram and density plot of junction amino acid length of clustered sequences from a sample
#'
#' @param clustered_seqs Clustered sequences from a sample
#'
#' @return junction amino acid length distribution from a sample
#' @export
len.sample.plot <- function(clustered_seqs) {
  nchar_junc <- nchar(clustered_seqs[, "junction_aa"])
  junc <- as.data.frame(nchar_junc)

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
      plot.title = element_text(lineheight = .8, size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14, colour = "black"),
      axis.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    )
}

#' Plot: Density ridges of junction amino acid length between clustered and unclustered sequence
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param clustered_seqs The clustered sequences
#' @param unclustered_seqs The unclustered sequences
#'
#' @return junction amino acid length distribution between the two groups
#' @export
len.clustered.plot <- function(clustered_seqs, unclustered_seqs) {
  df <- data.frame(
    group = c(
      rep("Clustered", nrow(clustered_seqs)),
      rep("Unclustered", nrow(unclustered_seqs))
    ),
    value = c(
      nchar(clustered_seqs$junction_aa),
      nchar(unclustered_seqs$junction_aa)
    )
  )

  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  print(paste0("p-value: ", res$p.format))
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(df, aes(x = value, y = group, color = group, fill = group)) +
    ggridges::geom_density_ridges(aes(color = group, fill = group),
      quantile_lines = TRUE, linewidth = 1.5, alpha = 0.75,
      vline_linetype = "dashed", quantiles = 1.5, bandwidth = 1
    ) +
    labs(x = "Length", y = "") +
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
    scale_fill_manual(values = c("#CCCC33", "#99FF99")) +
    scale_color_manual(values = c("#999900", "#66CC66")) +
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

#' Plot: Density ridges of junction amino acid length between two groups
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param group1_all_clustered_seqs All the clustered sequences in group1
#' @param group1_label Label of group1
#' @param group2_all_clustered_seqs All the clustered sequences in group2
#' @param group2_label Label of group2
#'
#' @return junction amino acid length distribution between the two groups
#' @export
len.group.plot <- function(group1_all_clustered_seqs, group1_label,
                           group2_all_clustered_seqs, group2_label) {
  df <- data.frame(
    group = c(
      rep(group1_label, nrow(group1_all_clustered_seqs)),
      rep(group2_label, nrow(group2_all_clustered_seqs))
    ),
    value = c(
      nchar(group1_all_clustered_seqs$junction_aa),
      nchar(group2_all_clustered_seqs$junction_aa)
    )
  )

  res <- ggpubr::compare_means(formula = value ~ group, df) %>% filter(p.signif != "ns")
  print(paste0("p-value: ", res$p.format))
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(df, aes(x = value, y = group, color = group, fill = group)) +
    ggridges::geom_density_ridges(aes(color = group, fill = group),
      quantile_lines = TRUE, linewidth = 1.5, alpha = 0.75,
      vline_linetype = "dashed", quantiles = 1.5, bandwidth = 1
    ) +
    labs(x = "Length", y = "") +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
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
  id <- bcr_clusters[[index]]$clonotype_index
  id <- id[loc]
  nn <- length(id)
  names <- c()
  for (i in 1:nn) {
    if (i < 10) {
      tmp.name <- paste("0", i, sep = "")
    } else {
      tmp.name <- i
    }
    names <- c(names, tmp.name)
  }

  return(names)
}

#' Plot: Visualization of multiple sequence alignment (MSA) of junction sequences within a cluster
#'
#' Perform MSA on the junction of all sequences in the clonal family for visualization
#'
#' @param bcr_clusters clonal families inferred by fastBCR
#' @param index Index of cluster
#' @param type 'DNA' for deoxyribonucleic acid or 'AA' for amino acid
#' @param raw_data The raw data which the clonal families inferred from
#' It is needed if you want to plot the MSA of 'DNA' sequences.
#' fastBCR will retrieve all the DNA sequences, which can be multiple sequences due to the degeneracy of codons, that correspond to the amino acid sequence of each clonotype from the raw data
#'
#' @return Visualization of MSA. From top to bottom are title('Vgene_Jgene_Length'), seqlogo, ggmsa, and msabar
#' @export
msa.plot <- function(bcr_clusters, index, type = c("AA", "DNA"), raw_data = NA) {
  if (type == "AA") {
    seqs_aa0 <- bcr_clusters[[index]]$junction_aa
    sink("default matrix.txt")
    MSAalign_aa <- msa::msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW", 
                            gapOpening = 20)
    sink()
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::AAStringSet(seqs_aa)
  } else if (type == "DNA") {
    raw_indices <- as.numeric(unique(unlist(strsplit(paste(bcr_clusters[[index]]$raw_indices, collapse = ','), ','))))
    seqs_aa0 <- raw_data$junction[raw_indices]
    sink("default matrix.txt")
    MSAalign_aa <- msa::msa(Biostrings::DNAStringSet(seqs_aa0), "ClustalW")
    sink()
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::DNAStringSet(seqs_aa)
  }

  seqs_aa_no <- gsub("-", "", seqs_aa)
  loc <- sort_msa(seqs_aa_no, seqs_aa0)
  SEQs_aa@ranges@NAMES <- name(bcr_clusters, index, loc)
  v_call <- unique(bcr_clusters[[index]]$v_call)
  j_call <- unique(bcr_clusters[[index]]$j_call)
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

#' Plot: Visualization of sequence logo of junction sequences within a cluster
#'
#' Perform sequence logo on the junction of all sequences in the clonal family for visualization
#'
#' @param bcr_clusters Clonal families inferred by fastBCR
#' @param index Index of cluster
#' @param type 'DNA' or 'AA'
#' @param raw_data The raw data which the clonal families inferred from
#' It is needed if you want to plot the sequence logo of 'DNA' sequences.
#' fastBCR will retrieve all the DNA sequences, which can be multiple sequences due to the degeneracy of codons, that correspond to the amino acid sequence of each clonotype from the raw data
#'
#' @return Sequence logo
#' @export
seqlogo.plot <- function(bcr_clusters, index, type = c("AA", "DNA"), raw_data = NA) {
  if (type == "AA") {
    seqs_aa0 <- bcr_clusters[[index]]$junction_aa
    sink("default matrix.txt")
    MSAalign_aa <- msa::msa(Biostrings::AAStringSet(seqs_aa0), "ClustalW", 
                            gapOpening = 20)
    sink()
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::AAStringSet(seqs_aa)
  } else if (type == "DNA") {
    raw_indices <- as.numeric(unique(unlist(strsplit(paste(bcr_clusters[[index]]$raw_indices, collapse = ','), ','))))
    seqs_aa0 <- raw_data$junction[raw_indices]
    sink("default matrix.txt")
    MSAalign_aa <- msa::msa(Biostrings::DNAStringSet(seqs_aa0), "ClustalW")
    sink()
    seqs_aa <- as.character(attributes(MSAalign_aa)$unmasked)
    SEQs_aa <- Biostrings::DNAStringSet(seqs_aa)
  }

  seqs_aa_no <- gsub("-", "", seqs_aa)
  loc <- sort_msa(seqs_aa_no, seqs_aa0)
  SEQs_aa@ranges@NAMES <- name(bcr_clusters, index, loc)
  v_call <- unique(bcr_clusters[[index]]$v_call)
  j_call <- unique(bcr_clusters[[index]]$j_call)
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

#' Function: Reconstructing a B cell lineage tree with minimum spanning tree and genotype abundances using ClonalTree
#'
#' The junction of BCR sequences within a cluster will be written as 'ClonalFamily_index.fasta' in 'ClonalTree/Examples/input' folder
#'
#' @param bcr_clusters Clonal families inferred by fastBCR
#' @param index Index of cluster
#' @param python_path The absolute path of the Python interpreter
#' @param raw_data The raw data which the clonal families inferred from
#' fastBCR will retrieve all the DNA sequences, which can be multiple sequences due to the degeneracy of codons, that correspond to the amino acid sequences of each clonotype from the raw data
#'
#' @return ClonalTree returns two files in the 'ClonalTree/Examples/output' folder
#' ClonalFamily_index.nk: the reconstructed BCR lineage tree in newick format
#' ClonalFamily_index.nk.csv: a table in csv format, containing the parent relationship and cost
#'
#' @export
clonal.tree.generation <- function(bcr_clusters, index, raw_data, python_path) {
  raw_indices <- as.numeric(unique(unlist(strsplit(paste(bcr_clusters[[index]]$raw_indices, collapse = ','), ','))))
  all_ids <- raw_data$sequence_id[raw_indices]
  all_junctions <- raw_data$junction[raw_indices]
  clonal_tree_df <- data.frame(sequence_id = all_ids, 
                               junction = all_junctions)
  all_junctions <- factor(all_junctions, levels = unique(all_junctions))
  clonal_tree_df <- clonal_tree_df[!duplicated(all_junctions), ]
  
  junction_seqs <- unique(clonal_tree_df$junction)
  MSAalign_dna <- msa::msa(Biostrings::DNAStringSet(junction_seqs), "ClustalW")
  align_junction_seqs <- as.character(attributes(MSAalign_dna)$unmasked)
  splt <- strsplit(align_junction_seqs, "")
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
  con.loc <- which(align_junction_seqs %in% consensus_DNA)
  naive_df <- data.frame(sequence_id = "naive", 
                         junction = consensus_DNA)
  if (length(con.loc) != 0) {
    clonal_tree_df <- clonal_tree_df[-con.loc, ]
  } 
  clonal_tree_df <- rbind(naive_df, clonal_tree_df)
  
  slash <- ifelse(grepl("\\\\", python_path), "\\", "/")
  file_name <- paste0("ClonalFamily_", index)
  fasta <- paste0(">", clonal_tree_df$sequence_id, "\n", clonal_tree_df$junction)
  write.table(fasta,
              file = paste0("ClonalTree", slash, "Examples", slash, "input", slash, file_name, ".fasta"),
              row.names = F,
              col.names = F,
              quote = F
  )
  
  system(paste0(python_path, " ClonalTree/src/clonalTree.py -i ClonalTree/Examples/input/", file_name, ".fasta -o ClonalTree/Examples/output/", file_name, ".abRT.nk -a 0 -r 1 -t 0"))
}

#' Plot: B cell lineage tree reconstructed by ClonalTree
#'
#' Orange points represent nodes and blue points represent tips
#' The x-axis shows the absolute genetic distance
#'
#' @param nk_path Path to the ClonalFamily_index.nk (the reconstructed BCR lineage tree in newick format)
#'
#' @return the reconstructed BCR lineage tree
#' @export
clonal.tree.plot <- function(nk_path) {
  tree <- ape::read.tree(nk_path)
  ggtree::ggtree(tree, size = 1, col = "deepskyblue3", options(ignore.negative.edge = TRUE)) +
    # ggtree::geom_nodelab(size = 4, color = "purple4", hjust = 1.2, vjust = -0.3) +
    ggtree::geom_nodepoint(size = 3, color = "orange2") +
    # ggtree::geom_tiplab(size = 4, color = "purple3", hjust = -0.3) +
    ggtree::geom_tippoint(size = 2, color = "deepskyblue3") +
    ggtree::theme_tree2() +
    theme(axis.text = element_text(face = "bold", size = 12, colour = "black"))
}

SHM.clu <- function(bcr_clusters, raw_data) {
  n <- length(bcr_clusters)
  ratio.lis <- rep(0, n)
  for (i in 1:n) {
    raw_indices <- as.numeric(unique(unlist(strsplit(paste(bcr_clusters[[i]]$raw_indices, collapse = ','), ','))))
    seqs <- raw_data$junction[raw_indices]
    seq.l <- max(nchar(seqs))
    nn <- length(seqs)
    clu.l <- seq.l * nn
    pair <- as.data.frame(t(combn(seqs, 2)))
    dist <- length(which(stringdist::stringdist(pair$V1, pair$V2, method = "lv") == 1))
    ratio.lis[i] <- dist / clu.l
  }

  return(ratio.lis)
}

SHM.iso <- function(bcr_clusters, raw_data) {
  n <- length(bcr_clusters)
  isotypes <- c("IGHD", "IGHM", "IGHA", "IGHG")
  all.dist <- rep(0, 4)
  all.len <- rep(0, 4)
  for (i in 1:n) {
    bcr_cluster <- bcr_clusters[[i]]
    iso <- substr(bcr_cluster$c_call, 1, 4)
    raw_indices <- as.numeric(unique(unlist(strsplit(paste(bcr_clusters[[index]]$raw_indices, collapse = ','), ','))))
    seqs <- raw_data$junction[raw_indices]

    for (is in 1:4) {
      tmp.loc <- which(iso == isotypes[is])
      if (length(tmp.loc) < 2) {
        next
      }
      iso.seqs <- seqs[tmp.loc]
      seq.l <- max(nchar(iso.seqs))
      nn <- length(iso.seqs)
      clu.l <- seq.l * nn
      pair <- as.data.frame(t(combn(iso.seqs, 2)))
      dist <- length(which(stringdist::stringdist(pair$V1, pair$V2, method = "lv") == 1))
      all.dist[is] <- all.dist[is] + dist
      all.len[is] <- all.len[is] + clu.l
    }
  }
  all.ratio <- all.dist / all.len

  return(all.ratio)
}

SHM.df <- function(clusters_list, raw_data_list, group_label) {
  average_SHM <- c()
  name <- c()

  pb <- utils::txtProgressBar(min = 0, max = length(clusters_list), style = 3)
  for (i in seq_along(clusters_list)) {
    clusters <- clusters_list[[i]]
    raw_data <- raw_data_list[[i]]
    SHM_clu <- SHM.clu(clusters, raw_data)
    SHM <- mean(SHM_clu)
    var_name <- names(clusters_list[i])
    name <- c(name, var_name)
    average_SHM <- c(average_SHM, SHM)
    setTxtProgressBar(pb, i)
  }

  SHM_df <- data.frame(
    group = rep(group_label, length(clusters_list)),
    SHM = average_SHM
  )
  rownames(SHM_df) <- name

  return(SHM_df)
}

SHM.iso.df <- function(clusters_list, raw_data_list, group_label) {
  SHM_iso <- c()
  names <- c()
  isotypes <- c("IGHD", "IGHM", "IGHA", "IGHG")

  pb <- utils::txtProgressBar(min = 0, max = length(clusters_list), style = 3)
  for (i in seq_along(clusters_list)) {
    clusters <- clusters_list[[i]]
    raw_data <- raw_data_list[[i]]
    SHM <- SHM.iso(clusters, raw_data)
    SHM_iso <- c(SHM_iso, SHM)
    var_name <- names(clusters_list[i])
    names <- c(names, paste(var_name, isotypes, sep = "_"))
    setTxtProgressBar(pb, i)
  }

  SHM_iso_df <- data.frame(
    group = rep(group_label, 4 * length(clusters_list)),
    SHM = SHM_iso,
    isotypes = rep(isotypes, length(clusters_list))
  )
  rownames(SHM_iso_df) <- names

  return(SHM_iso_df)
}

#' Function: Calculate the average SHM ratio of all clusters and the SHM ratios of clustered sequences in four isotypes in each sample
#'
#' @param clusters_list1 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group1
#' @param raw_data_list1 A list where each element is the raw data named after its filename in group1
#' @param group1_label Label of group1
#' @param clusters_list2 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group2
#' @param raw_data_list2 A list where each element is the raw data named after its filename in group2
#' @param group2_label Label of group2
#'
#' @return Dataframe containing the average SHM ratio and SHM ratios in four isotypes
#' @export
SHM.calculation <- function(clusters_list1, raw_data_list1, group1_label,
                            clusters_list2, raw_data_list2, group2_label) {
  print("group1 processing:")
  df1 <- SHM.df(clusters_list1, raw_data_list1, group1_label)
  print("group2 processing:")
  df2 <- SHM.df(clusters_list2, raw_data_list2, group2_label)
  SHM_df <- rbind(df1, df2)

  return(SHM_df)
}

#' Function: Calculate the the SHM ratios of clustered sequences in four isotypes in each sample
#'
#' @param clusters_list1 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group1
#' @param raw_data_list1 A list where each element is the raw data named after its filename in group1
#' @param group1_label Label of group1
#' @param clusters_list2 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group2
#' @param raw_data_list2 A list where each element is the raw data named after its filename in group2
#' @param group2_label Label of group2
#'
#' @return Dataframe containing the SHM ratios in four isotypes
#' @export
SHM.iso.calculation <- function(clusters_list1, raw_data_list1, group1_label,
                                clusters_list2, raw_data_list2, group2_label) {
  print("group1 processing:")
  df1 <- SHM.iso.df(clusters_list1, raw_data_list1, group1_label)
  print("group2 processing:")
  df2 <- SHM.iso.df(clusters_list2, raw_data_list2, group2_label)
  SHM_iso_df <- rbind(df1, df2)

  return(SHM_iso_df)
}

#' Plot: Boxplot of the SHM ratios between the two groups
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param SHM_df Dataframe containing the average SHM ratio of all clusters and the SHM ratios of clustered sequences in four isotypes in each sample
#'
#' @return Boxplot showing the SHM ratios between the two groups
#' @export
SHM.plot <- function(SHM_df) {
  res <- ggpubr::compare_means(formula = SHM ~ group, SHM_df) %>% filter(p.signif != "ns")
  print(paste0("p-value: ", res$p.format))
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = SHM_df, aes(x = group, y = SHM, fill = group, color = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6)) +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T
    ) +
    coord_cartesian(clip = "off")
}

#' Plot: Boxplot of SHM ratios of clustered sequences in different isotypes
#'
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param SHM_df Dataframe containing the average SHM ratio of all clusters and the SHM ratios of clustered sequences in four isotypes in each sample
#'
#' @return Boxplot showing the SHM ratios of clustered sequences in four isotypes ('IGHD', 'IGHM', 'IGHA', and 'IGHG') between two groups
#' @export
SHM.iso.plot <- function(SHM_iso_df) {
  SHM_iso_df$group_iso <- paste(SHM_iso_df$group, SHM_iso_df$isotypes, sep = "_")
  res <- ggpubr::compare_means(formula = SHM ~ group_iso, SHM_iso_df) %>% filter(p.signif != "ns")
  print(paste0(res$group1, "~", res$group2, " p-value: ", res$p.format))
  my_comparisons <- list()
  ii <- 1
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      tmp1 <- unlist(strsplit(res$group1[i], "_"))[2]
      tmp2 <- unlist(strsplit(res$group2[i], "_"))[2]
      if (tmp1 == tmp2) {
        my_comparisons[[ii]] <- c(res$group1[i], res$group2[i])
        ii <- ii + 1
      }
    }
  }

  ggplot(data = SHM_iso_df, aes(
    x = factor(group_iso, levels = c(
      paste(unique(SHM_iso_df$group), "IGHD", sep = "_"),
      paste(unique(SHM_iso_df$group), "IGHM", sep = "_"),
      paste(unique(SHM_iso_df$group), "IGHA", sep = "_"),
      paste(unique(SHM_iso_df$group), "IGHG", sep = "_")
    )), y = SHM,
    fill = group,
    color = group
  )) +
    stat_boxplot(aes(color = group), geom = "errorbar", width = 0.1, linewidth = 0.5, position = position_dodge(0.6)) +
    geom_boxplot(aes(fill = group), position = position_dodge(0.6), linewidth = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = group), position = position_jitterdodge(), shape = 21, size = 1, alpha = 0.9) +
    labs(x = "", y = "SHM ratio") +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    scale_x_discrete(labels = c(
      "IGHM", "IGHM",
      "IGHD", "IGHD",
      "IGHA", "IGHA",
      "IGHG", "IGHG"
    )) +
    coord_cartesian(clip = "off") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = TRUE,
      step.increase = 0.1
    )
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

#' Plot: Network of CSR within a cluster
#'
#' @param bcr_clusters Clonal families inferred by fastBCR
#' @param index Index of cluster
#'
#' @return Visualization of isotype co-occurrence within a BCR cluster.
#' Circle size represents the number of sequences carrying a given isotype.
#' Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses.
#' The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
#' Matrix of values of connected edges between clustered sequences in different isotypes is printed.
#' @export
CSR.cluster.plot <- function(bcr_clusters, index) {
  cluster_name <- deparse(substitute(bcr_clusters))
  title_name <- paste0(cluster_name, ": index = ", index)

  tmp <- bcr_clusters[[index]]
  bcr_clusters.count <- nrow(tmp)
  tmp.cw <- matrix(0, 1, ncol = 11)
  colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
  tmp.cw[, "Cluster"] <- index
  tmp.is <- tmp$c_call
  tmp.is <- fastBCR:::Iso.fst(tmp.is)
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
  print(network.matrix)
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
    guides(size = "none") +
    ggtitle(title_name)
}

#' Plot: Network of CSR within a sample
#'
#' @param bcr_clusters Clonal families inferred by fastBCR
#'
#' @return Visualization of isotype co-occurrence within BCR clusters. Circle size represents the number of clusters carrying a given isotype.
#' Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses.
#' The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among clusters.
#' Matrix of values of connected edges between clustered sequences in different isotypes is printed.
#' @export
CSR.sample.plot <- function(bcr_clusters) {
  title_name <- deparse(substitute(bcr_clusters))

  bcr.cluster.isotypes <- NULL
  for (i in seq_along(bcr_clusters)) {
    tmp <- bcr_clusters[[i]]
    tmp.cw <- matrix(0, 1, ncol = 11)
    colnames(tmp.cw) <- c("Cluster", "IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2", "Unidentified")
    tmp.cw[, "Cluster"] <- i
    tmp.is <- tmp$c_call
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
  print(network.matrix)
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
    guides(size = "none") +
    ggtitle(title_name)
}

NAb.seq.query <- function(pro_data, AbDab, method, maxDist) {
  join <- inner_join(pro_data, AbDab, by = c("v_call", "j_call"), relationship = "many-to-many")
  colnames(join)[c((ncol(pro_data) + 1):ncol(join))] <- paste0("NAb_", colnames(join)[c((ncol(pro_data) + 1):ncol(join))])
  colnames(join)[which(colnames(join) == "CDRH3.x")] <- "CDRH3"
  colnames(join)[which(colnames(join) == "NAb_CDRH3.y")] <- "NAb_CDRH3"

  if (!is.na(method) && !is.na(maxDist)) {
    result <- join %>%
      filter(stringdist::stringdist(CDRH3, NAb_CDRH3, method = method) <= maxDist)
  } else {
    # Exact matching when neither method nor maxDist is given
    result <- join %>%
      filter(stringr::str_detect(stringr::fixed(CDRH3), stringr::fixed(NAb_CDRH3)))
  }

  NAb_n <- length(unique(result$clonotype_index))

  return(NAb_n)
}

NAb.ratio <- function(pro_data, bcr_clusters, AbDab, method, maxDist, species) {
  # AbDab
  AbDab_splits <- t(sapply(AbDab$Heavy.V.Gene, NAb.Species.splits))
  AbDab$v_call <- AbDab_splits[, 1]
  AbDab$v_species <- AbDab_splits[, 2]
  AbDab_splits <- t(sapply(AbDab$Heavy.J.Gene, NAb.Species.splits))
  AbDab$j_call <- AbDab_splits[, 1]
  AbDab$j_species <- AbDab_splits[, 2]
  AbDab <- AbDab[tolower(AbDab$v_species) == tolower(species), ]
  AbDab <- AbDab[tolower(AbDab$j_species) == tolower(species), ]
  AbDab <- AbDab[complete.cases(AbDab[, c("v_call", "j_call")]), ]
  AbDab <- select(AbDab, v_call, j_call, CDRH3) %>%
    filter(nchar(CDRH3) >= 7)

  # pro_data
  pro_data$CDRH3 <- substr(pro_data$junction_aa, 2, nchar(pro_data$junction_aa) - 1)
  pro_data <- select(pro_data, clonotype_index, v_call, j_call, CDRH3)
  pro_NAb_n <- NAb.seq.query(pro_data, AbDab, method, maxDist)

  if (pro_NAb_n == 0) {
    NAb_ratio <- 0
  } else {
    clustered_index <- unique(unlist(sapply(bcr_clusters, function(x) x$clonotype_index)))
    clustered_seqs <- pro_data[clustered_index, ]
    clu_NAb_n <- NAb.seq.query(clustered_seqs, AbDab, method, maxDist)
    NAb_ratio <- clu_NAb_n / pro_NAb_n
  }

  return(NAb_ratio)
}

NAb.ratio.df <- function(pro_data_list, clusters_list, group_label, AbDab, method, maxDist, species) {
  sample_n <- length(pro_data_list)
  NAb_ratios <- rep(0, sample_n)
  name <- rep(0, sample_n)
  pb <- utils::txtProgressBar(min = 0, max = sample_n, style = 3)
  for (i in 1:sample_n) {
    pro_data <- pro_data_list[[i]]
    clusters <- clusters_list[[i]]
    NAb_ratios[i] <- NAb.ratio(
      pro_data, clusters,
      AbDab, method, maxDist, species
    )
    name[i] <- names(clusters_list[i])
    setTxtProgressBar(pb, i)
  }
  NAb_ratio_df <- data.frame(
    group = rep(group_label, length(clusters_list)),
    NAb_ratio = NAb_ratios
  )
  rownames(NAb_ratio_df) <- name

  return(NAb_ratio_df)
}

#' Function: Calculate the NAb ratios between the two groups.
#'
#' NAb ratio is established as an indicator of the proportional prevalence of neutralizing antibody sequences within expanded clonotypes in each sample
#' It is defined as the fraction of the number of NAb sequences within clonal families to the total number of NAb sequences present in each sample
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param pro_data_list1 A list where each element is the preprocessed data named after its filename in group1
#' @param clusters_list1 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group1
#' @param group1_label Label of group1
#' @param pro_data_list2 A list where each element is the preprocessed data named after its filename in group2
#' @param clusters_list2 A list where each element of the list is the clonal families inferred by fastBCR for each sample in group2
#' @param group2_label Label of group2
#' @param NAb_vjcdr3 IGHV gene, IGHJ gene and CDRH3 of all NAb sequences
#'
#' @return Dataframe containing the NAb ratios between the two groups.
#' @export
NAb.ratio.calculation <- function(pro_data_list1, clusters_list1, group1_label,
                                  pro_data_list2, clusters_list2, group2_label,
                                  AbDab, method = NA, maxDist = NA, species = "Human") {
  print("group1 processing:")
  df1 <- NAb.ratio.df(
    pro_data_list1, clusters_list1, group1_label,
    AbDab, method, maxDist, species
  )
  print("group2 processing:")
  df2 <- NAb.ratio.df(
    pro_data_list2, clusters_list2, group2_label,
    AbDab, method, maxDist, species
  )
  NAb_ratio_df <- rbind(df1, df2)

  return(NAb_ratio_df)
}

#' Plot: Boxplot of the NAb ratios between the two groups.
#'
#' NAb ratio is established as an indicator of the proportional prevalence of neutralizing antibody sequences within expanded clonotypes in each sample
#' It is defined as the fraction of the number of NAb sequences within clonal families to the total number of NAb sequences present in each sample
#' Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test
#'
#' @param NAb_ratio_df Dataframe containing the NAb ratios between the two groups.
#'
#' @return Boxplot showing the NAb ratios between the two groups.
#' @export
NAb.ratio.plot <- function(NAb_ratio_df) {
  res <- ggpubr::compare_means(formula = NAb_ratio ~ group, NAb_ratio_df) %>% filter(p.signif != "ns")
  print(paste0("p-value: ", res$p.format))
  my_comparisons <- list()
  if (nrow(res) != 0) {
    for (i in 1:nrow(res)) {
      my_comparisons[[i]] <- c(res$group1[i], res$group2[i])
    }
  }
  ggplot(data = NAb_ratio_df, aes(x = group, y = NAb_ratio, fill = group, color = group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5, position = position_dodge(0.6)) +
    geom_boxplot(position = position_dodge(0.6), size = 0.5, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(fill = group), position = position_jitter(0.1), shape = 21, size = 2, alpha = 0.9) +
    xlab("") +
    ylab("NAb ratio") +
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
    scale_fill_manual(values = c("#FF3333", "#0099CC")) +
    scale_color_manual(values = c("#CC0000", "#336699")) +
    coord_cartesian(clip = "off") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      hide.ns = T,
      step.increase = 0.1
    ) +
    scale_y_continuous(limits = c(0, 1))
}

# Split the columns and remove the parentheses
NAb.Species.splits <- function(x) {
  # Determine whether it contains 'ND' first, and return NA if it does
  if (grepl("ND", x)) {
    return(c(NA, NA))
  }

  splits <- strsplit(x, " \\(")[[1]]
  if (length(splits) > 1) {
    v_call <- splits[[1]]
    species <- gsub("\\)", "", splits[[2]])
  } else {
    return(c(NA, NA))
  }

  return(c(v_call, species))
}

#' Function: Query the corresponding sequence from the public antibody database
#'
#' @param bcr_clusters Clonal families inferred by fastBCR
#' @param AbDab The public antibody database
#' @param method CDRH3 matching method. It can be 'NA' for perfect match, 'hamming' for hamming distance or 'lv' for Levenshtein distance. Defaults to 'NA'.
#' @param maxDist Maximum distance allowed for matching when the argument 'method' is 'hamming' or 'lv'. Defaults to 'NA'.
#' @param species can be 'Mouse' or 'Human'. Defaults to 'Human'.
#'
#' @return Result of the query
#' @export
NAb.query <- function(bcr_clusters, AbDab, method = NA, maxDist = NA, species = "Human") {
  # Record the associated cluster id and merge all the clusters data frames
  clusters_df <- bind_rows(bcr_clusters, .id = "cluster_id")
  clusters_df$CDRH3 <- substr(clusters_df$junction_aa, 2, nchar(clusters_df$junction_aa) - 1)

  # Data preprocessing:
  AbDab_splits <- t(sapply(AbDab$Heavy.V.Gene, NAb.Species.splits))
  # Split the v_call and species
  AbDab$v_call <- AbDab_splits[, 1]
  AbDab$v_species <- AbDab_splits[, 2]
  AbDab_splits <- t(sapply(AbDab$Heavy.J.Gene, NAb.Species.splits))
  # Split the j_call and species
  AbDab$j_call <- AbDab_splits[, 1]
  AbDab$j_species <- AbDab_splits[, 2]
  AbDab <- AbDab[tolower(AbDab$v_species) == tolower(species), ]
  AbDab <- AbDab[tolower(AbDab$j_species) == tolower(species), ]
  AbDab <- AbDab[complete.cases(AbDab[, c("v_call", "j_call")]), ]

  # Intra-join clusters and dfs to match exactly according to v_call and j_call
  join <- inner_join(clusters_df, AbDab, by = c("v_call", "j_call"), relationship = "many-to-many")
  colnames(join)[c((ncol(clusters_df) + 1):ncol(join))] <- paste0("NAb_", colnames(join)[c((ncol(clusters_df) + 1):ncol(join))])
  colnames(join)[which(colnames(join) == "CDRH3.x")] <- "CDRH3"
  colnames(join)[which(colnames(join) == "NAb_CDRH3.y")] <- "NAb_CDRH3"

  if (!is.na(method) && !is.na(maxDist)) {
    result <- join %>%
      filter(stringdist::stringdist(CDRH3, NAb_CDRH3, method = method) <= maxDist)
  } else {
    # Exact matching when neither method nor maxDist is given
    result <- join %>%
      filter(stringr::str_detect(CDRH3, stringr::fixed(NAb_CDRH3)) & stringr::str_detect(NAb_CDRH3, stringr::fixed(CDRH3)))
  }

  if (nrow(result) > 0) {
    print(paste("The query result has", nrow(result), "rows."))
    result <- select(
      result, cluster_id, clonotype_index, v_call, j_call, junction_aa,
      NAb_Name, NAb_Heavy.V.Gene, NAb_Heavy.J.Gene, NAb_CDRH3
    )
  } else {
    print("The query result is empty.")
  }

  return(result)
}
