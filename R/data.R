#' simulated single cell multi-omics data.
#'
#' A dataset containing paired single-cell RNA-seq and ATAC-seq data as well as the true labels of cell groups
#'
#' @format A list contains two field data:
#' \describe{
#'   \item{data}{Simulated single-cell RNA-seq and ATAC-seq data}
#'   \item{labels}{The true labels of cell groups}
#'   ...
#' }
#' @source \url{https://bitbucket.org/ConesaLab/mosim/src/master/}
"pairedData_simulation"



#' single cell multi-omics data of A549 cells.
#'
#' A549 dataset containing paired single-cell RNA-seq and ATAC-seq data as well as the collected time information of each cell
#'
#' @format A list contains two field data:
#' \describe{
#'   \item{data}{Single-cell RNA-seq and ATAC-seq data}
#'   \item{labels}{The time information of cell groups}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089}
"pairedData_A549"
