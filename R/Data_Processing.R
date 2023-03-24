#' @title Function: The first step of fastBCR to process raw data.
#'
#' @description The input of the function needs to meet the AIRR standard format (containing at least 'sequence_id', 'v_call', 'j_call', and 'junction_aa' information). Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved. Sequences with the same 'v_call', 'j_call' and 'junction_aa' are considered identical and deduplicated.
#'
#' @param raw_data AIRR format data. Inference of clonal familes requires following columns to be present in the raw_data: "sequence_id", "v_call", "j_call", "junction_aa". "junction" or "c_call" are optional if you want to plot "evolutionary tree" or need isotypes related analyis (SHM/CSR).
#'
#' @return Processed data as input for BCR_Clustering
#' @export
#'
#' @examples
#' data("test_data")
#' input <- data.pro(test_data)
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

  return(pro_data)
}
