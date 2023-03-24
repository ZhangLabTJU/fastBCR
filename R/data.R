#' Data: As a sample dataset for Data_Processing.
#'
#' A real experimental BCR sequencing repertoire with AIRR format.
#'
#' @examples
#' data("test_data")
#' input = data.pro(test_data)
"test_data"

#' Data: As a sample dataset for BCR_clustering.
#'
#' Processed data after Data_Processing (using 'data.pro(test_data)').
#'
#' @examples
#' data("input")
#' bcr_clusters = BCR.cluster(input)
"input"

#' Data: As a sample dataset for Evolutionary_Analysis.
#'
#' Inferred clonal familes by fastBCR (using 'BCR.cluster(input)').
#'
#' @examples
#' data("bcr_clusters")
"bcr_clusters"

#' Data: As a sample data for plot.SHM.sample (sample1).
#'
#' Example1 for SHM ratios of four isotypes (using 'SHM.sample(bcr_clusters_1)').
#'
#' @examples
#' data("SHM_ratio_1")
"SHM_ratio_1"

#' Data: As a sample data for plot.SHM.sample (sample2).
#'
#' Example2 for SHM ratios of four isotypes (using 'SHM.sample(bcr_clusters_2)').
#'
#' @examples
#' data("SHM_ratio_2")
"SHM_ratio_2"

#' Data: As a sample data for plot.SHM.sample (sample3).
#'
#' Example3 for SHM ratios of four isotypes (using 'SHM.sample(bcr_clusters_3)').
#'
#' @examples
#' data("SHM_ratio_3")
"SHM_ratio_3"

#' Data: As a sample data for plot.SHM.sample (sample4).
#'
#' Example4 for SHM ratios of four isotypes (using 'SHM.sample(bcr_clusters_4)').
#'
#' @examples
#' data("SHM_ratio_4")
"SHM_ratio_4"

#' Data: As a sample data for plot.SHM.sample (sample5).
#'
#' Example5 for SHM ratios of four isotypes (using 'SHM.sample(bcr_clusters_5)').
#'
#' @examples
#' data("SHM_ratio_5")
"SHM_ratio_5"

#' ighv_hum_df
#'
#' human germline IgH-V (heavy chain v gene segments). When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 119 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"ighv_hum_df"


#' ighd_hum_df
#'
#' human germline IgH-D (heavy chain d-gene segments). When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 37 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"ighd_hum_df"


#' ighj_hum_df
#'
#' human germline IgH-V (heavy chain v-gene segments). When multiple alleles were present, the
#' first one was included. These names and sequences can be changed by
#' customized by changing this dataframe. Additionally, repeating elements
#' can give certain germline gene elements a larger probability of being used
#' during repertoire evolution.
#'
#' @format A data frame with 6 rows and 2 variables:
#' \describe{
#'   \item{gene}{The gene name}
#'   \item{seq}{The corresponding sequence}
#' }
#' @source IMGT
"ighj_hum_df"
