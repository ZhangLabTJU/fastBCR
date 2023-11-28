#' Data: As a sample data for Simulation
#'
"germline_data"

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
