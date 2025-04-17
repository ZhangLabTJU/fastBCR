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
  raw_data_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(file_list), style = 3)
  for (i in seq_along(file_list)) {
    file_name <- basename(file_list[i])
    var_name <- unlist(strsplit(file_name, "\\."))[1]
    compressed <- sum(sapply(archive_format, function(x) length(grep(x, file_name))))
    if (compressed) {
      if (storage_format == "csv") {
        raw_data_list[[var_name]] <- read.table(archive::archive_read(file_list[i]), header = TRUE, sep = ",")
      } else if (storage_format == "tsv") {
        raw_data_list[[var_name]] <- read.table(archive::archive_read(file_list[i]), header = TRUE, sep = "\t")
      } else if (storage_format == "Rdata") {
        load(archive::archive_read(file_list[i]))
        loaded_objects <- load(archive::archive_read(file_list[i]))
        raw_data_list[[var_name]] <- get(loaded_objects)
      }
    } else {
      if (storage_format == "csv") {
        raw_data_list[[var_name]] <- read.table(file_list[i], header = TRUE, sep = ",")
      } else if (storage_format == "tsv") {
        raw_data_list[[var_name]] <- read.table(file_list[i], header = TRUE, sep = "\t")
      } else if (storage_format == "Rdata") {
        load(file_list[i])
        loaded_objects <- load(file_list[i])
        raw_data_list[[var_name]] <- get(loaded_objects)
      }
    }
    setTxtProgressBar(pb, i)
  }
  
  cat(paste0(
    paste0("\nYou have loaded ", length(raw_data_list), " samples\n"),
    paste(paste0("Sample '", names(raw_data_list), "' contains ", sapply(raw_data_list, function(x) nrow(x)), " sequences"), collapse = "\n")
  ))
  
  return(raw_data_list)
}

data.productive <- function(raw_data, productive_only){
  # filter sequences with appropriate length
  productive_data <- dplyr::filter(raw_data, !is.na(v_call) & !is.na(j_call) & !is.na(junction_aa)) %>% # filter sequences without 'v_call', 'j_call' or 'junction_aa' entry
    dplyr::filter(nchar(junction_aa) >= 9 & nchar(junction_aa) <= 26) # filter sequences with too short or too long junction length
  
  # productive only
  if (productive_only) {
    prod_col <- grep("^productive$", names(productive_data), ignore.case = TRUE, value = TRUE)
    if (length(prod_col) > 0) {
      productive_data <- productive_data[tolower(as.character(productive_data[[prod_col]])) != "false", ]
      productive_data <- productive_data[tolower(as.character(productive_data[[prod_col]])) != "f", ]
    }
  }
  junction_aa <- productive_data$junction_aa
  star <- grep("\\*", junction_aa)
  underline <- grep("_", junction_aa)
  X <- grep("X", junction_aa)
  rm <- unique(c(star, underline, X))
  if (length(rm) != 0) {
    productive_data <- productive_data[-rm, ]
  } else {
    productive_data <- productive_data
  }
  
  # V/J gene without allele
  v_call <- strsplit(productive_data$v_call, "\\*")
  v_call <- sapply(v_call, function(x) x[1])
  j_call <- strsplit(productive_data$j_call, "\\*")
  j_call <- sapply(j_call, function(x) x[1])
  productive_data$v_call <- v_call
  productive_data$j_call <- j_call
  
  return(productive_data)
}

#' @title Function: The first step of fastBCR to process raw data.
#'
#' @description The input data needs to contain essential columns including "v_call" (V gene with or without allele), "j_call" (J gene with or without allele) and "junction_aa" (amino acid translation of the junction).
#' Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved. You can also filter clonotype sequence with low frequency to reduce sequencing errors to some extent.
#' Sequences with the same "v_call", "j_call" and "junction_aa" are considered to be the same clonotype and are merged into one row in processed data.
#' The column "clonotype_index" is the index of each clonotype.
#' The column "clonotype_count" is the count of each clonotype.
#' The column "clone_count" is the sum of the counts (calculated based on "count_col_name" parameter) of the sequences belonging to each clonotype.
#' The column "clone_fre" is the frequency version of "clone_count".
#' The column "raw_index" is the index of each clonotype representative sequence in raw data.
#' The column "raw_indices" is the indices of all sequences with each clonotype in raw data.
#'
#' @param raw_data Raw data for preprocessed.
#' Inference of clonal familes requires following columns to be present in the raw_data: "v_call" (V gene with or without allele), "j_call" (J gene with or without allele) and "junction_aa" (amino acid translation of the junction).
#' The "junction" (junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons) column is needed for phylogenetic tree construction and SHM-related analysis
#' The "c_call" (constant region gene with or without allele) column is needed for isotype-related analysis.
#' @param productive_only Whether to filter out unproductive sequences based on the "productive" information in the raw data.
#'  Defaults to "FALSE".
#' @param count_col_name The column name for the count of each sequence.
#'  Defaults to "NA" which means the original count of the sequence is not taken into account.
#' @param count_filter_thre The threshold of "clone_count" under which the clonotype sequence will be filtered.
#'  Defaults to "NA" which means no clonotype sequence should be filtered.
#'
#' @return Processed data as input for clonal family inference
#' @export
data.pro <- function(raw_data, productive_only = FALSE, count_col_name = NA, count_filter_thre = NA) {
  raw_data$raw_index <- c(1:nrow(raw_data))
  # check if the raw_data meet the requirements on column names
  required_col_names <- c("v_call", "j_call", "junction_aa")
  if (all(required_col_names %in% colnames(raw_data)) == FALSE) {
    stop("The required column names ('v_call', 'j_call', 'junction_aa') must be present in the raw_data. Please check the format of the input data or change the column name to meet the requirements.")
  }
  
  # get productive data with appropriate length
  productive_data <- data.productive(raw_data, productive_only)
  productive_data <- data.frame(productive_data)
  
  # clonotype (v-j-junction_aa) deduplicate
  v_call <- productive_data$v_call
  j_call <- productive_data$j_call
  junction_aa <- productive_data$junction_aa
  if (!is.na(count_col_name)) {
    count <- productive_data[, count_col_name]
    count[is.na(count)] <- 1
  } else {
    count <- rep(1, nrow(productive_data))
  }
  v_j_junction_aa <- paste(v_call, j_call, junction_aa)
  v_j_junction_aa <- factor(v_j_junction_aa, levels = unique(v_j_junction_aa))
  pro_data <- productive_data[!duplicated(v_j_junction_aa), ]
  
  if(length(pro_data) == 0){
    return(NA)
  }else{
    # clonotype backtrack
    index_match <- data.frame(
      clonotype_numeric = as.numeric(v_j_junction_aa),
      clonotype_index = 1:length(v_j_junction_aa)
    )
    index_match <- I(split(index_match[, 2], list(index_match$clonotype_numeric), drop = TRUE))
    clonotype_n <- length(unique(v_j_junction_aa))
    clone_count <- rep(0, clonotype_n)
    clonotype_count <- rep(0, clonotype_n)
    max_freq_index <- rep(0, clonotype_n)
    raw_indices <- rep(0, clonotype_n)
    
    for (i in 1:clonotype_n) {
      tmp.index <- index_match[[i]]
      clone_count[i] <- sum(count[tmp.index])
      clonotype_count[i] <- length(tmp.index)
      max_freq_index[i] <- tmp.index[which.max(count[tmp.index])]
      raw_indices[i] <- paste(productive_data[tmp.index, 'raw_index'], collapse = ',')
    }
    
    pro_data <- data.frame(
      clonotype_index = c(1:clonotype_n), clonotype_count, clone_count,
      clone_fre = clone_count / sum(clone_count),
      productive_data[max_freq_index, ], raw_indices
    )
    if (!is.na(count_filter_thre)) {
      count <- pro_data$clone_count
      pro_data <- dplyr::filter(pro_data, count >= count_filter_thre)
    }
    if(length(pro_data) == 0){
      return(NA)
    }
    rownames(pro_data) <- 1:nrow(pro_data)
    
    return(pro_data)
  }
}

#' Function: Preprocessing of raw data to meet the input requirements for clonal family inference
#'
#' @param raw_data_list A list where each element is the raw data named after its filename
#' @param productive_only Whether to filter out unproductive sequences based on the "productive" information in the raw data.
#'  Defaults to "FALSE".
#' @param count_col_name The column name for the count of each sequence.
#'  Defaults to "NA" which means the original count of the sequence is not taken into account.
#' @param count_filter_thre The threshold of "clone_count" under which the clonotype sequence will be filtered.
#'  Defaults to "NA" which means no clonotype sequence should be filtered.
#'
#' @return A list where each element is the processed data named after its filename
#' @export
data.preprocess <- function(raw_data_list, productive_only = FALSE, count_col_name = NA, count_filter_thre = NA) {
  pro_data_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(raw_data_list), style = 3)
  for (i in seq_along(raw_data_list)) {
    raw_data <- raw_data_list[[i]]
    var_name <- names(raw_data_list[i])
    processed_data <- data.pro(raw_data, productive_only, count_col_name, count_filter_thre)
    pro_data_list[[var_name]] <- processed_data
    setTxtProgressBar(pb, i)
  }
  
  cat(paste0(
    paste0("\nYou have preprocessed ", length(pro_data_list), " samples\n"),
    paste(paste0("Sample '", names(pro_data_list), "' contains ", sapply(pro_data_list, function(x) nrow(x)), " unique clonotypes"), collapse = "\n")
  ))
  
  return(pro_data_list)
}