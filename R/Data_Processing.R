#' @title Function: Preprocessing of raw data to meet fastBCR requirements for input data
#'
#' @description The input of the function needs to meet the AIRR standard format (containing at least 'sequence_id', 'v_call', 'j_call', and 'junction_aa' information). Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved. Sequences with the same 'v_call', 'j_call' and 'junction_aa' are considered identical and deduplicated.
#'
#' @param raw_data AIRR format data. Inference of clonal familes requires following columns to be present in the raw_data: 'sequence_id', 'v_call', 'j_call', 'junction_aa'. 'junction' or 'c_call' are optional if you want to plot 'evolutionary tree' or need isotypes related analyis (SHM/CSR).
#'
#' @return Processed data as input for clonal family inference
#' @export
data.pro <- function(raw_data) {
  # sequence_id unique
  seq_ids <- raw_data$sequence_id
  raw_data <- raw_data[!duplicated(seq_ids), ]

  # junction length
  loc <- which(nchar(raw_data$junction_aa) >= 9 & nchar(raw_data$junction_aa) <= 26)
  raw_data <- raw_data[loc, ]
  not.triple <- which(nchar(raw_data$junction) %% 3 != 0)
  if (length(not.triple) != 0) {
    raw_data <- raw_data[-not.triple, ]
  }

  # productive only
  junction_aa <- raw_data$junction_aa
  star <- grep("\\*", junction_aa)
  underline <- grep("_", junction_aa)
  X <- grep("X", junction_aa)
  v_call_na <- which(is.na(raw_data$v_call))
  j_call_na <- which(is.na(raw_data$j_call))
  rm <- c(star, underline, X, v_call_na, j_call_na)
  if (length(rm) != 0) {
    productive_data <- raw_data[-rm, ]
  } else {
    productive_data <- raw_data
  }

  # v-j-junction unique
  v_call <- strsplit(productive_data$v_call, "\\*")
  v_call <- sapply(v_call, function(x) x[1])
  j_call <- strsplit(productive_data$j_call, "\\*")
  j_call <- sapply(j_call, function(x) x[1])
  productive_data$v_call <- v_call
  productive_data$j_call <- j_call
  junction_aa <- productive_data$junction_aa
  v_j_junction <- paste(v_call, j_call, junction_aa)
  v_j_junction <- factor(v_j_junction, levels = unique(v_j_junction))
  pro_data <- productive_data[!duplicated(v_j_junction), ]
  tab <- table(v_j_junction)
  pro_data$clonotype_ids <- paste("clonotype", c(1:nrow(pro_data)), sep = "_")
  pro_data$clonotype_number <- as.numeric(tab)
  pro_data$clonotype_frequency <- as.numeric(tab) / nrow(productive_data)

  return(pro_data)
}
