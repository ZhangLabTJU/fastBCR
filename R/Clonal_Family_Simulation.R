data.pro <- function(raw_data) {
  # junction length
  loc <- which(nchar(raw_data$junction_aa) >= 9 & nchar(raw_data$junction_aa) <= 26)
  pro_data <- raw_data[loc, ]
  not.triple <- which(nchar(pro_data$junction) %% 3 != 0)
  if (length(not.triple) != 0) {
    pro_data <- pro_data[-not.triple, ]
  }

  # productive only
  junction_aa <- pro_data$junction_aa
  star <- grep("\\*", junction_aa)
  underline <- grep("_", junction_aa)
  X <- grep("X", junction_aa)
  unpro <- union(star, underline)
  rm <- union(unpro, X)
  if (length(rm) != 0) {
    pro_data <- pro_data[-rm, ]
  }

  # v-j-junction unique
  rm <- union(which(is.na(pro_data$v_call)), which(is.na(pro_data$j_call)))
  if (length(rm) != 0) {
    pro_data <- pro_data[-rm, ]
  }
  v_call <- strsplit(pro_data$v_call, "\\*")
  v_call <- sapply(v_call, function(x) x[1])
  j_call <- strsplit(pro_data$j_call, "\\*")
  j_call <- sapply(j_call, function(x) x[1])
  pro_data$v_call <- v_call
  pro_data$j_call <- j_call
  junction_aa <- pro_data$junction_aa
  vjcdr3 <- paste(v_call, j_call, junction_aa)
  pro_data <- pro_data[!duplicated(vjcdr3), ]
  seq_ids <- pro_data$sequence_id
  pro_data <- pro_data[!duplicated(seq_ids), ]
  return(pro_data)
}

# data.frame to fasta
df2fas <- function(data, filename) {
  ids <- data$ids
  seqs <- data$seqs
  fasta <- paste0(">", ids, "\n", seqs)
  write.table(fasta,
    file = filename,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# Insertion
ins <- function(dna, loc) {
  left_dna <- substring(dna, 1, loc)
  len <- nchar(dna)
  right_dna <- substring(dna, loc + 1, len)
  ran_ins_dna <- sample(x = c("A", "T", "G", "C"), 1, replace = TRUE, c(rep(1 / 4, 4)))
  ins_dna <- paste(left_dna, ran_ins_dna, right_dna, sep = "")
  return(ins_dna)
}

# Deletion
del <- function(dna, loc) {
  left_dna <- substring(dna, 1, loc - 1)
  len <- nchar(dna)
  right_dna <- substring(dna, loc + 1, len)
  del_dna <- paste(left_dna, right_dna, sep = "")
  return(del_dna)
}

# Substitution
sub <- function(dna, loc) {
  holding_char <- substr(dna, loc, loc)
  if (holding_char == "A") {
    substr(dna, loc, loc) <- sample(x = c("T", "G", "C"), 1, replace = TRUE, c(15, 70, 15))
  } else if (holding_char == "T") {
    substr(dna, loc, loc) <- sample(x = c("A", "G", "C"), 1, replace = TRUE, c(15, 15, 70))
  } else if (holding_char == "G") {
    substr(dna, loc, loc) <- sample(x = c("A", "T", "C"), 1, replace = TRUE, c(70, 15, 15))
  } else if (holding_char == "C") {
    substr(dna, loc, loc) <- sample(x = c("A", "G", "T"), 1, replace = TRUE, c(15, 15, 70))
  }
  return(dna)
}

#' Function: Generation of a fasta file consisting of germline sequences.
#'
#' @description Randomly generate germline sequences based in V(D)J recombination and convert them to a fasta file.
#'
#' @param gemline_num The number of simulated germline sequences.
#' @param filename Path to save fasta file.
#'
#' @return A saved fasta file consisting of randomly generated germline sequences that need to be annotated subsequently.
#' @export
#'
#' @examples
#' germline2fas(100, filename = "Simulation/Germline.fasta")
germline2fas <- function(gemline_num, filename) {
  data("ighv_hum_df")
  data("ighd_hum_df")
  data("ighj_hum_df")
  indV <- sample(x = seq_len(nrow(ighv_hum_df)), gemline_num, replace = TRUE)
  indD <- sample(x = seq_len(nrow(ighd_hum_df)), gemline_num, replace = TRUE)
  indJ <- sample(x = seq_len(nrow(ighj_hum_df)), gemline_num, replace = TRUE)
  vseq <- toupper(as.character(ighv_hum_df[[2]][indV]))
  dseq <- toupper(as.character(ighd_hum_df[[2]][indD]))
  jseq <- toupper(as.character(ighj_hum_df[[2]][indJ]))

  # D-J rearranged
  djseq <- c()
  for (i in 1:gemline_num) {
    tmp.dseq <- dseq[i]
    dseq_length <- nchar(tmp.dseq)
    dseq_left <- stringr::str_sub(tmp.dseq, 1, round(dseq_length / 2))
    dseq_right <- stringr::str_sub(tmp.dseq, round(dseq_length / 2) + 1, dseq_length)
    tmp.jseq <- jseq[i]
    jseq_length <- nchar(tmp.jseq)
    jseq_left <- stringr::str_sub(tmp.jseq, 1, 15)
    jseq_right <- stringr::str_sub(tmp.jseq, 16, jseq_length)
    junction <- paste(dseq_right, jseq_left, sep = "")

    # 1-3 Indels randomly occur
    indel_num <- sample(1:3, 1)
    for (j in 1:indel_num) {
      indel_mode <- sample(1:2, 1)
      loc <- sample(1:nchar(junction), 1)
      if (indel_mode == 1) {
        junction <- ins(junction, loc)
      } else {
        junction <- del(junction, loc)
      }
    }
    tmp.djseq <- paste(dseq_left, junction, jseq_right, sep = "")
    djseq <- c(djseq, tmp.djseq)
  }

  # V-D-J rearranged
  vdjseq <- c()
  for (i in 1:gemline_num) {
    tmp.vseq <- vseq[i]
    vseq_length <- nchar(tmp.vseq)
    vseq_left <- stringr::str_sub(tmp.vseq, 1, (vseq_length - 15))
    vseq_right <- stringr::str_sub(tmp.vseq, (vseq_length - 14), vseq_length)
    tmp.djseq <- djseq[i]
    djseq_length <- nchar(tmp.djseq)
    djseq_left <- stringr::str_sub(tmp.djseq, 1, 15)
    djseq_right <- stringr::str_sub(tmp.djseq, 16, djseq_length)
    junction <- paste(vseq_right, djseq_left, sep = "")

    # 1-3 Indels randomly occur
    indel_num <- sample(1:3, 1)
    for (j in 1:indel_num) {
      indel_mode <- sample(1:2, 1)
      loc <- sample(1:nchar(junction), 1)
      if (indel_mode == 1) {
        junction <- ins(junction, loc)
      } else {
        junction <- del(junction, loc)
      }
    }
    tmp.vdjseq <- paste(vseq_left, junction, djseq_right, sep = "")
    vdjseq <- c(vdjseq, tmp.vdjseq)
  }

  ids <- c(1:gemline_num)
  seqs <- vdjseq
  fasta.df <- data.frame(ids, seqs)
  df2fas(fasta.df, filename)
}

each_round <- function(input, mut, start, end) {
  mode <- c(rep("ins", 4), rep("del", 5), rep("sub", 500)) # 0.8:1:100
  mut_n <- length(input)
  not_mut <- 1 - mut
  mut_dna <- c()
  for (kk in 1:mut_n) {
    cdr3_dna <- input[kk]
    act_n <- sample(3:5, 1)
    act_dna <- rep(cdr3_dna, act_n)
    for (ii in 1:act_n) {
      tmp_dna <- act_dna[ii]
      mut_loc <- c()
      mut_form <- c()
      for (cc in start:end) {
        CDR_mut <- sample(x = c(0, 1), 1, replace = TRUE, c(not_mut, mut))
        if (CDR_mut == 1) {
          form <- mode[sample(1:509, 1)]
          mut_loc <- c(mut_loc, cc)
          mut_form <- c(mut_form, form)
        }
      }
      if (length(mut_loc) == 0) {
        mut_dna <- c(mut_dna, tmp_dna)
      } else {
        for (ll in seq_along(length(mut_loc))) {
          tmp.form <- mut_form[ll]
          tmp.loc <- mut_loc[ll]
          if (tmp.form == "ins") {
            tmp_dna <- ins(tmp_dna, tmp.loc)
          } else if (tmp.form == "del") {
            tmp_dna <- del(tmp_dna, tmp.loc)
          } else {
            tmp_dna <- sub(tmp_dna, tmp.loc)
          }
        }
        mut_dna <- c(mut_dna, tmp_dna)
      }
    }
  }
  drop_off <- 0.4
  drop_num <- round(drop_off * length(mut_dna))
  mut_dna <- mut_dna[-sample(1:length(mut_dna), drop_num, replace = FALSE)]
  return(mut_dna)
}

#' Function: Generation of a fasta file consisting of sequences from distinct clonal familes.
#'
#' @param germline_data Annotated germline sequences.
#' @param CF_n The number of simulated clonal families.
#' @param mut_ratio The point mutation rate for SHM during simulated affinity maturation process whose selectable range is (0,1).
#' @param filename Path to save fasta file.
#'
#' @return A saved fasta file consisting of sequences from distinct clonal familes that need to be annotated subsequently.
#' @export
#'
#' @examples
#' data("germline_data") # read.table('Simulation/Germline_igblast_db-pass_parse-select.tsv', header = T, sep = "\t")
#' CF2fas(germline_data, CF_n = 10, mut_ratio = 0.001, filename = 'Simulation/10_0.001.fasta')

CF2fas <- function(germline_data, CF_n, mut_ratio, filename) {
  germline_data <- data.pro(germline_data)
  data <- germline_data[sample(seq_len(nrow(germline_data)), CF_n, replace = FALSE), ]
  final_dna <- NULL
  mode <- c(rep("ins", 4), rep("del", 5), rep("sub", 500)) # 0.8:1:100
  cdr3_start <- data$v_sequence_end - 15
  cdr3_end <- data$j_sequence_start + 15
  not_mut_ratio <- 1 - mut_ratio

  for (mm in 1:CF_n) {
    # first round
    ori_dna <- data[mm, "sequence"]
    start <- cdr3_start[mm]
    end <- cdr3_end[mm]
    act_n <- 5
    act_dna <- rep(ori_dna, act_n)
    mut_dna1 <- c()
    for (ii in 1:act_n) {
      tmp_dna <- act_dna[ii]
      mut_loc <- c()
      mut_form <- c()
      for (cc in start:end) {
        CDR_mut <- sample(x = c(0, 1), 1, replace = TRUE, c(not_mut_ratio, mut_ratio))
        if (CDR_mut == 1) {
          form <- mode[sample(1:509, 1)]
          mut_loc <- c(mut_loc, cc)
          mut_form <- c(mut_form, form)
        }
      }
      if (length(mut_loc) == 0) {
        mut_dna1 <- c(mut_dna1, tmp_dna)
      } else {
        for (ll in seq_along(length(mut_loc))) {
          tmp.form <- mut_form[ll]
          tmp.loc <- mut_loc[ll]
          if (tmp.form == "ins") {
            tmp_dna <- ins(tmp_dna, tmp.loc)
          } else if (tmp.form == "del") {
            tmp_dna <- del(tmp_dna, tmp.loc)
          } else {
            tmp_dna <- sub(tmp_dna, tmp.loc)
          }
        }
        mut_dna1 <- c(mut_dna1, tmp_dna)
      }
    }
    drop_off <- 0.4
    drop_num <- round(drop_off * length(mut_dna1))
    mut_dna1 <- mut_dna1[-sample(1:length(mut_dna1), drop_num, replace = FALSE)]

    mut_dna2 <- each_round(mut_dna1, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])
    mut_dna3 <- each_round(mut_dna2, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])
    mut_dna4 <- each_round(mut_dna3, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])
    mut_dna5 <- each_round(mut_dna4, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])
    mut_dna6 <- each_round(mut_dna5, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])
    mut_dna7 <- each_round(mut_dna6, mut = mut_ratio, start = cdr3_start[mm], end = cdr3_end[mm])

    clu <- c(ori_dna, mut_dna1, mut_dna2, mut_dna3, mut_dna4, mut_dna5, mut_dna6, mut_dna7)
    final_dna[[mm]] <- clu
  }

  # final_dna to fasta
  ids <- c()
  seqs <- c()
  for (i in 1:CF_n) {
    tmp.seqs <- final_dna[[i]]
    seqs <- c(seqs, tmp.seqs)
    for (j in 1:length(tmp.seqs)) {
      tmp.ids <- paste(i, "_", j, sep = "")
      ids <- c(ids, tmp.ids)
    }
  }
  fasta.df <- data.frame(ids, seqs)
  df2fas(fasta.df, filename)
}
