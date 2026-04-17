#' atacInferCnv: CNV inference from single-cell chromatin accessibility data
#'
#' \code{atacInferCnv} provides tools for genome-wide inference, analysis,
#' and visualization of copy number variation (CNV) from single cell
#' and single nucleus ATAC-seq data. The package extends the use of
#' \code{infercnv} for chromatin accessibility profiles and includes utilities
#' for preparing input matrices, aggregating genomic regions into bins,
#' constructing meta-cells, running CNV inference, and visualizing inferred
#' CNV profiles.
#'
#' The package is intended for the analysis of tumor single-cell chromatin
#' accessibility datasets and supports identification of
#' large-scale chromosomal alterations, and subclonal CNV structure.
#'
#' @section Key functions:
#' \describe{
#'   \item{\code{\link{prepareAtacInferCnvInput}}}{Prepare input matrices and
#'   annotations for CNV inference from scATAC-seq data.}
#'   \item{\code{\link{runAtacInferCnv}}}{Run CNV calling.}
#'   \item{\code{\link{plotCnvBlocks}}}{Visualize inferred CNV block profiles.}
#' }
#'
#' @return The package-level help page is used for documentation only and does
#'   not return a value.
#' @name atacInferCnv
#' @aliases atacInferCnv-package
#' @keywords package
#' @importFrom utils write.csv read.delim write.table read.table
#' @importFrom config get
#' @importFrom stats cutree
#' @importFrom utils capture.output
#' @import graphics
#' @import grDevices
#' @import ggplot2
#' @import Seurat
#' @import Signac
#' @import GenomicRanges
#' @import SingleCellExperiment
#' @import infercnv
#' @import stringr
#' @useDynLib atacInferCnv, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom SummarizedExperiment assay colData
NULL
