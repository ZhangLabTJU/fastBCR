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
#' The column "clonotype_count" is the count of each clonotype.
#' The column "clone_count" is the sum of the counts (calculated based on "count_col_name" parameter) of the sequences belonging to each clonotype.
#' The column "clone_fre" is the frequency version of "clone_count".
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
    stop("The required column names ('sequence_id', 'v_call', 'j_call', 'junction_aa') must be present in the raw_data. Please XXX")
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

  # V/J gene without allele
  v_call <- strsplit(productive_data$v_call, "\\*")
  v_call <- sapply(v_call, function(x) x[1])
  j_call <- strsplit(productive_data$j_call, "\\*")
  j_call <- sapply(j_call, function(x) x[1])
  productive_data$v_call <- v_call
  productive_data$j_call <- j_call

  # clonotype (v-j-junction_aa) deduplicate
  junction_aa <- productive_data$junction_aa
  if(!is.na(count_col_name)){
    count <- productive_data[,count_col_name]
  }else{
    count <- rep(1, nrow(productive_data))
  }
  v_j_junction_aa <- paste(v_call, j_call, junction_aa)
  v_j_junction_aa <- factor(v_j_junction_aa, levels = unique(v_j_junction_aa))
  pro_data <- productive_data[!duplicated(v_j_junction_aa), ]

  # clonotype backtrack
  index_match <- data.frame(clonotype_index = as.numeric(v_j_junction_aa),
                            orign_index = 1:length(v_j_junction_aa))
  index_match <- I(split(index_match[,2], list(index_match$clonotype_index), drop = TRUE))
  clonotype_n <- length(unique(v_j_junction_aa))
  clone_count = rep(0, clonotype_n)
  clonotype_count = rep(0, clonotype_n)
  max_freq_index  = rep(0, clonotype_n)

  # pb <- utils::txtProgressBar(min = 0, max = clonotype_n, style = 3)
  for(i in 1:clonotype_n){
    tmp.index <- index_match[[i]]
    tmp.clone_count <- sum(count[tmp.index])
    tmp.clonotype_count <- length(tmp.index)
    clone_count[i] = tmp.clone_count
    clonotype_count[i] = tmp.clonotype_count
    max_freq_index[i] = tmp.index[which.max(count[tmp.index])]
    # setTxtProgressBar(pb, i)
  }

  pro_data <- data.frame(clonotype_index = c(1:clonotype_n), clonotype_count, clone_count,
                         clone_fre = clone_count/sum(clone_count), orign_index = max_freq_index,
                         productive_data[max_freq_index,], index_match)
  rownames(pro_data) = 1:nrow(pro_data)

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
set <- c('a','a','b','c','d','a','c','b','a','c')
set <- data.frame(index= 1:length(set), num=set)
index.match = data.frame(orign.index = set$index, new.index = as.numeric(factor(set$num)))
index.match$orign.index[which(index.match$new.index == 1)]
