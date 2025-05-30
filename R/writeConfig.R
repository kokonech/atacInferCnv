#' Write InferCNV configuration
#'
#' @param resDir Result directory path
#' @param sId Result name
#' @param ctrlGrp Name for the reference control cell type
#' @param binSize Size of the bin e.g. 100000 for 100Kbp
#' @param meta True if use meta cells, default false
#'
#' @return NULL
#' @export
#'
writeConfig <- function(resDir, sId , ctrlGrp = "Normal", binSize=NULL,meta=F ) {

  if (meta) {
    inputFile = paste0(sId,"_raw_matrix_metacells.txt.gz")
    annFile = paste0(sId,"_annotation_metacells.tsv")
    cutOff = 0.25
  } else {
    inputFile = paste0(sId,"_raw_counts.txt.gz")
    annFile = paste0(sId,"_cnv_ann.txt")
    cutOff = 0.1
  }

  customRef = paste0(sId, "_cnv_ref.txt")

  if (!is.null(binSize)) {
    customRef = gsub(".txt",sprintf(".binsize_%d.txt", binSize), customRef)
  }


  numClusters = 3

  config_content <- paste0("
  default:
    resName: '", paste0(sId,"_infercnv"),"'
    countsFile: '", inputFile, "'
    annFile: '", annFile, "'
    customRef: '", customRef, "'
    refGroups: '", ctrlGrp, "'
    cutOff: ", cutOff, "
    numClusters: ", numClusters, "
  ")

  # Write to config.yml
  writeLines(config_content, paste0(resDir,"infercnv_config.yml"))
}
