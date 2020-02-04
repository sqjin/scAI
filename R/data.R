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
#' @source \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1932-8}
"data_simulation"

#' A simulation dataset consists of five imbalanced cell clusters with five clusters in scRNA-seq data and three clusters in scATAC-seq data
#'
#' @format A list contains two field data:
#' \describe{
#'   \item{data}{Simulated single-cell RNA-seq and ATAC-seq data}
#'   \item{labels}{The true labels of cell groups}
#'   ...
#' }
#' @source \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1932-8}
"data_simulation8"


#' single cell multi-omics data of kidney cells.
#'
#' kidney dataset containing paired single-cell RNA-seq and ATAC-seq data as well as the label information of each cell from the original paper
#'
#' @format A list contains two field data:
#' \describe{
#'   \item{data}{Single-cell RNA-seq and ATAC-seq data}
#'   \item{labels}{The cell type information of cell groups}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089}
"data_kidney"


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
"data_A549"


#' single cell multi-omics data of mESC cells.
#'
#' mESC dataset containing paired single-cell RNA-seq and single-cell methylation data as well as the culture condition information of each cell
#'
#' @format A list contains two field data:
#' \describe{
#'   \item{data}{Single-cell RNA-seq and single-cell methylation data}
#'   \item{labels}{The culture condition information of cells}
#'   ...
#' }
#' @source \url{https://github.com/bioFAM/MOFA/blob/master/vignettes/MOFA_example_scMT.Rmd}
"data_mESC"

