#' Function: Load datasets into a list
#'
#' The compressed files in "7zip", "cab", "cpio", "iso9660", "lha", "mtree", "shar", "rar", "raw", "tar", "tar.gz", "xar", "zip", "warc" format can also be read in
#'
#' @param folder_path Path to the data folder
#' @param storage_format The storage format of data. It can be in "csv", "tsv", or "Rdata" format.
#'
#' @return A list where each element is the raw data named after its filename
#' @export
data.load <- function(folder_path, storage_format) {
  info <- file.info(folder_path)
  if (info$isdir) {
    file_list <- list.files(folder_path, full.names = TRUE)
  } else {
    file_list <- folder_path
  }
  archive_format <- c("7zip", "cab", "cpio", "iso9660", "lha", "mtree", "shar", "rar", "raw", "tar", "tar.gz", "xar", "zip", "warc")
  data_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(file_list), style = 3)
  for (i in seq_along(file_list)) {
    file_name <- basename(file_list[i])
    var_name <- unlist(strsplit(file_name, "\\."))[1]
    compressed <- sum(sapply(archive_format, function(x) length(grep(x, file_name))))
    if (compressed) {
      if (storage_format == "csv") {
        data_list[[var_name]] <- read.table(archive::archive_read(file_list[i]), header = TRUE, sep = ",")
      } else if (storage_format == "tsv") {
        data_list[[var_name]] <- read.table(archive::archive_read(file_list[i]), header = TRUE, sep = "\t")
      } else if (storage_format == "Rdata") {
        load(archive::archive_read(file_list[i]))
        data_list[[var_name]] <- get(var_name)
      }
    } else {
      if (storage_format == "csv") {
        data_list[[var_name]] <- read.table(file_list[i], header = TRUE, sep = ",")
      } else if (storage_format == "tsv") {
        data_list[[var_name]] <- read.table(file_list[i], header = TRUE, sep = "\t")
      } else if (storage_format == "Rdata") {
        load(file_list[i])
        data_list[[var_name]] <- get(var_name)
      }
    }
    setTxtProgressBar(pb, i)
  }

  return(data_list)
}

#' @title Function: The first step of fastBCR to process raw data.
#'
#' @description The input of the function needs to meet the AIRR standard format (containing at least "sequence_id", "v_call", "j_call", and "junction_aa" information).
#' Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved.
#' Sequences with the same "v_call", "j_call" and "junction_aa" are considered to be the same clonotype and are merged into one row in processed data.
#' In each row of processed data, "sequence_id" is the "sequence_id" of sequences with the same clonotype separated by ";".
#' "clonotype_count" and "clonotype_fre" is the count and frequency of the clonotype calculated based on "count_col_name" parameter.
#' If "junction" ("c_call") is contained in the raw data, it is replaced by the most frequently occurring "junction" ("c_call") in the clonotype.
#'
#' @param raw_data AIRR format data.
#' Inference of clonal familes requires following columns to be present in the raw_data: "sequence_id", "v_call", "j_call", "junction_aa".
#' "junction" or "c_call" are optional if you want to plot the evolutionary tree or need isotypes related (SHM/CSR) analysis.
#' @param count_col_name The column name for the count of each sequence.
#'  It can be "consensus_count", "duplicate_count" or "umi_count" according to your needs.
#'  Defaults to "NA" which means the original count of the sequence is not taken into account.
#'
#' @return Processed data as input for clonal family inference
#' @export
data.pro <- function(raw_data, count_col_name = NA) {
  # check if the raw_data meet the requirements on column names
  required_col_names <- c("sequence_id", "v_call", "j_call", "junction_aa")
  if(all(required_col_names %in% colnames(raw_data)) == FALSE){
    stop("The required column names ('sequence_id', 'v_call', 'j_call', 'junction_aa') must be present in the raw_data")
  }

  # junction length filtering
  raw_data <- dplyr::filter(raw_data, nchar(junction_aa) >= 9 & nchar(junction_aa) <= 26)

  # productive only
  junction_aa <- raw_data$junction_aa
  star <- grep("\\*", junction_aa)
  underline <- grep("_", junction_aa)
  X <- grep("X", junction_aa)
  v_call_na <- which(is.na(raw_data$v_call))
  j_call_na <- which(is.na(raw_data$j_call))
  rm <- unique(c(star, underline, X, v_call_na, j_call_na))
  if (length(rm) != 0) {
    productive_data <- raw_data[-rm, ]
  } else {
    productive_data <- raw_data
  }

  # clonotype (v-j-junction_aa) deduplicate
  v_call <- strsplit(productive_data$v_call, "\\*")
  v_call <- sapply(v_call, function(x) x[1])
  j_call <- strsplit(productive_data$j_call, "\\*")
  j_call <- sapply(j_call, function(x) x[1])
  productive_data$v_call <- v_call
  productive_data$j_call <- j_call
  sequence_id <- productive_data$sequence_id
  junction_aa <- productive_data$junction_aa
  if(!is.na(count_col_name)){
    counts <- productive_data[,count_col_name]
  }
  if(any(colnames(productive_data) %in% "junction")) junction <- productive_data$junction
  if(any(colnames(productive_data) %in% "c_call")) c_call <- productive_data$c_call
  v_j_junction_aa <- paste(v_call, j_call, junction_aa)
  v_j_junction_aa <- factor(v_j_junction_aa, levels = unique(v_j_junction_aa))
  pro_data <- productive_data[!duplicated(v_j_junction_aa), ]

  # clonotype backtrack
  index_match <- data.frame(clonotype_index = as.numeric(v_j_junction_aa),
                            orign_index = 1:length(v_j_junction_aa))
  clonotype_n <- length(unique(v_j_junction_aa))
  clonotype_sequence_id <- c()
  clonotype_counts <- c()
  if(any(colnames(pro_data) %in% "junction")) clonotype_junction <- c()
  if(any(colnames(pro_data) %in% "c_call")) clonotype_c_call <- c()
  for(i in 1:clonotype_n){
    tmp.index <- index_match$orign_index[which(index_match$clonotype_index == i)]
    tmp.clonotype_sequence_id <- paste(sequence_id[tmp.index], collapse = ";")
    clonotype_sequence_id <- c(clonotype_sequence_id, tmp.clonotype_sequence_id)
    if(!is.na(count_col_name)){
      tmp.clonotype_counts <- sum(counts[tmp.index])
    }else{
      tmp.clonotype_counts <- length(tmp.index)
    }
    clonotype_counts <- c(clonotype_counts, tmp.clonotype_counts)
    if(any(colnames(pro_data) %in% "junction")){
      tmp.clonotype_junction <- names(sort(table(junction[tmp.index]),decreasing = T))[1]
      clonotype_junction <- c(clonotype_junction, tmp.clonotype_junction)
    }
    if(any(colnames(pro_data) %in% "c_call")){
      tmp.c_call <- names(sort(table(c_call[tmp.index]),decreasing = T))[1]
      clonotype_c_call <- c(clonotype_c_call, tmp.c_call)
    }
  }

  pro_data$clonotype_index <- c(1:clonotype_n)
  pro_data$sequence_id <- clonotype_sequence_id
  pro_data$clonotype_count <- clonotype_counts
  pro_data$clonotype_fre <- clonotype_counts/sum(clonotype_counts)
  if(any(colnames(pro_data) %in% "junction")) pro_data$junction <- clonotype_junction
  if(any(colnames(pro_data) %in% "c_call")) pro_data$c_call <- clonotype_c_call

  pro_data <- dplyr::select(pro_data, clonotype_index, sequence_id, v_call, j_call,
                     junction_aa, clonotype_count, clonotype_fre, everything())

  return(pro_data)
}

#' Function: Preprocessing of raw data to meet the input requirements for clonal family inference
#'
#' @param data_list A list where each element is the raw data named after its filename
#' @param count_col_name The column name for the count of each sequence.
#'  It can be "consensus_count", "duplicate_count" or "umi_count" according to your needs.
#'  Defaults to "NA" which means the original count of the sequence is not taken into account.
#'
#' @return A list where each element is the processed data named after its filename
#' @export
data.preprocess <- function(data_list, count_col_name = NA) {
  pro_data_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(data_list), style = 3)
  for (i in seq_along(data_list)) {
    data <- data_list[[i]]
    var_name <- names(data_list[i])
    processed_data <- data.pro(data, count_col_name)
    pro_data_list[[var_name]] <- processed_data
    setTxtProgressBar(pb, i)
  }

  return(pro_data_list)
}
