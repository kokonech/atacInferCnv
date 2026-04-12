
pa_message <- function(msg) {
  # source: https://ropensci.org/blog/2024/02/06/verbosity-control-packages/
  if (getOption("perfromGA.quiet", FALSE)) {
    return()
  }
  message(msg)
}


mergeWithControl <- function(mb, ctrlObj, targColumn) {
  # by default number of control cells <= 33% of number of tumor cells
  # TODO: make a param
  expNumCtrlCells <- round( 0.33 * ncol(mb))
  if (ncol(ctrlObj) > expNumCtrlCells) {
    pa_message(paste0("Adjust external control, decrease num cells to ",
                      expNumCtrlCells))
    ctrlObj <- ctrlObj[ , seq_len(expNumCtrlCells)]
  }
  #print(ctrlObj)
  #print(ctrlObj@assays$ATAC@ranges)
  #print(mb@assays$ATAC@ranges)
  # NOTE : not most optimal merging - regions of x are used as main

  mb.merged <- merge(mb, y=ctrlObj)

  mb.merged@meta.data[ , targColumn] <- c(mb@meta.data[,targColumn],
                                          rep("ExtControl",ncol(ctrlObj)) )
  mb.merged

}


prepareAnalysis <- function(mb, resDir, sId, annData, targColumn,
                            ctrlObj, performGA, paVerbose) {
  cIds <- intersect(rownames(annData), rownames(mb@meta.data))
  if (length(cIds) == 0) {
    stop("No overlapping cells IDs in annotation!")
  }
  pa_message("Adjust for custom annotation")
  mb <- subset(mb, cells= cIds)
  #mb <- mb[ ,cIds]
  #print(summary(rownames(mb@meta.data) == cIds))
  mb@meta.data <- cbind(mb@meta.data,
                        annData[rownames(mb@meta.data),targColumn])
  colnames(mb@meta.data)[ncol(mb@meta.data)] <- targColumn

  if (!is.null(ctrlObj)) {
    pa_message("Merge with external control object")
    mb <- mergeWithControl(mb, ctrlObj, targColumn )
  }


  if (performGA) {
    pa_message("Perform general analysis...")
    pa_message("Normalization...")

    mb <- RunTFIDF(mb, verbose = paVerbose)
    mb <- FindTopFeatures(mb, min.cutoff = 'q0', verbose = paVerbose)
    mb <- RunSVD(mb, verbose = paVerbose)

    pa_message("Dimensional reduction...")

    ndim <- 30 # default 30

    mb <- RunUMAP(object = mb, reduction = 'lsi', dims = 2:ndim,
                  verbose = paVerbose)
    mb <- FindNeighbors(object = mb, reduction = 'lsi', dims = 2:ndim,
                        verbose = paVerbose)
    mb <- FindClusters(object = mb, verbose = paVerbose, algorithm = 3)

    pdf(file.path(resDir,paste0(sId,"_UMAP.pdf")),width = 8, height = 6)
    # print REQUIRED for ggplot2 output
    print(DimPlot(object = mb, pt.size=1, label=TRUE))
    if (nchar(targColumn) > 0) {
      # REQUIRED for ggplot2 output
      print(DimPlot(mb, reduction = "umap",
                    label = TRUE,group.by = targColumn))
    }
    dev.off()
  }
  saveRDS(mb, file.path(resDir,paste0(sId,"_obj.RDS") ))

  mb
}


saveCnvInput <- function(mb,resDir, sId, targColumn) {

  # for InferCNV

  annTable <- mb@meta.data
  annTable2 <- annTable[,targColumn,drop=FALSE]

  #if (targColumn == "seurat_clusters" ) {
  #    annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
  #}
  write.table(annTable2, file.path(resDir,paste0(sId,"_cnv_ann.txt")) ,
              col.names = FALSE,sep="\t",quote=FALSE)

  gz1 <- gzfile(file.path(resDir,paste0(sId,"_raw_counts.txt.gz")), "w")
  rawCounts <- as.matrix(mb@assays$ATAC@counts)

  rawCounts <- rawCounts[ grep("chr",rownames(rawCounts)),]
  write.table(rawCounts,gz1,sep="\t",quote=FALSE)
  close(gz1)

  peakRegions <- GRanges(sub("-",":",rownames(rawCounts), fixed=TRUE))
  seqlevels(peakRegions) <- paste0("chr",c(seq_len(22),"X","Y")) # fix order
  peakRegions <- sort(peakRegions)


  peakDf <- data.frame(peakRegions)[,seq_len(3)]
  rownames(peakDf) <- paste0(  peakDf[,1],"-",peakDf[,2],"-",peakDf[,3] )
  peakDf$seqnames <- gsub("chr","", peakDf$seqnames)
  summary(rownames(rawCounts) %in% rownames(peakDf))
  write.table(peakDf, file.path(resDir,paste0(sId,"_cnv_ref.txt")),
              col.names = FALSE,sep="\t",quote=FALSE)
}

#' Prepare input for the CNV calling from scATAC-seq data
#'
#' @param dataPath Path to the input data in 10X format or ATAC counts matrix
#'  in txt format (gzipped)
#' @param annPath Path to annotation of the cells in tab-delimited format
#' @param resDir Path to the result directory
#' @param inObj Pre-computed Seurat/Signac object with required input data
#'  (alternative for dataPath)
#' @param sId Result name. Default: "Sample"
#' @param targColumn Name of the target column in annotation. Default:CellType
#' @param ctrlGrp Name for the reference control cell type. Could be several
#' names, separated by comma. Default: "Normal"
#' @param ctrlObj Seurat/Signac object to use as non-tumor control. Default:NULL
#' @param binSize Apply custom bin size to combine signals in windows for CNV
#' calling e.g. 500000 for 500 KBp. Default: NULL (not use this option)
#' @param chromLength Numeric vector of chromosome sizes, specific for genome.
#'  Default: NULL
#' @param metaCells Set TRUE to use meta cells (n=5 cells by default) or assign
#'  a number of cells. Default: FALSE
#' @param performGA Perform general analysis of scATAC-data (clustering, UMAP).
#'  Default: TRUE
#' @param verbose Detailed output, progress messages, and diagnostic
#' information. Default: TRUE
#' @return Invisibly returns NULL.
#' @examples
#' resPath = tempfile()
#' inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz",
#'                       package = "atacInferCnv")
#' sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann_n30.txt",
#'                      package = "atacInferCnv" )
#' prepareAtacInferCnvInput(inPath,sAnn,resPath, targColumn = "cnvBlock",
#'                           ctrlGrp = "Normal",performGA = FALSE)
#'
#' @export
#'
prepareAtacInferCnvInput <- function(dataPath = "",
                                     annPath = "",
                                     resDir = "", inObj = NULL, sId = "sample",
                                     targColumn = "CellType",ctrlGrp = "Normal",
                                     ctrlObj = NULL,
                                     binSize = NULL, chromLength = NULL,
                                     metaCells = FALSE,
                                     performGA = TRUE, verbose = TRUE) {



  # take care about messages
  options(perfromGA.quiet = !verbose)

  pa_message("Loading input...")
  if (nchar(resDir) == 0) {
    stop("Path to result is not provided!")
  }
  if (is.null(inObj)) {
    if (!(file.exists(dataPath))) {
      stop("Input data folder is not found:",dataPath)

    }
    if (dir.exists(dataPath)) {
      pa_message("Input is directory, assuming 10X data folder...")
      # scMulti-omics
      countsPath <- file.path(dataPath,"filtered_feature_bc_matrix.h5")
      fragpath <- file.path(dataPath, "atac_fragments.tsv.gz")
      # scAtac
      countsPath2 <- file.path(dataPath,"filtered_peak_bc_matrix.h5")
      if (!(file.exists(countsPath))) {
        if (file.exists(countsPath2)) {
          pa_message("10X scATAC format identified")
          countsVals <- Read10X_h5(countsPath2)
          fragpath <- file.path(dataPath, "fragments.tsv.gz")
        } else {
          stop("Input feature counts matrix is not found in path: ",dataPath,
             "\nExpected formats: filtered_feature_bc_matrix.h5
             or filtered_peak_bc_matrix.h5")
        }
      } else {
        pa_message("10X scMulti-omics format identified")
        countsVals <- Read10X_h5(countsPath)$Peaks
      }


      if (!(file.exists(fragpath))) {
        stop("Input fragments loci is not found:",fragpath)
      }

      # create ATAC assay and add it to the object
      chrom_assay  <- CreateChromatinAssay(
        counts = countsVals,
        sep = c(":", "-"),
        fragments = fragpath
      )
    } else {
      pa_message("Input is a file, assuming peaks matrix...")
      countsVals <- read.delim(dataPath,check.names = FALSE)
      chrom_assay  <- CreateChromatinAssay(
        counts = countsVals,
        sep = c(":", "-")
      )
    }

    mb <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "ATAC",
      project = sId
    )

    if (!(file.exists(annPath))) {
      stop("Annotation file is not found: ",annPath)

    }
    annData <- read.delim(annPath)

    if (! (targColumn %in% colnames(annData) ) ) {
      stop("Required annotation column is not available: ", targColumn)
    } else {
      pa_message(paste0("Using target annotation column: ",targColumn))
      #print(class(annData[, targColumn]))
      if (is.null(ctrlObj)) {
        annInfo <- summary(as.factor(annData[, targColumn]))
        pa_message(paste(capture.output(annInfo)))
        ctrlStatus <- as.numeric(
          str_split_1(ctrlGrp,pattern = ",") %in% names(annInfo)
          )
        if (any(ctrlStatus == 0) ) {
          stop("Non-tumor control group is not found in annotation: ", ctrlGrp)
        }
      } else {
        if (!inherits(ctrlObj, "Seurat")) {
          stop("Non-tumor external control input is not Seurat object!")
        }
        pa_message("Using external control (assigned as ExtControl)")
        #print(ctrlObj)
        ctrlGrp <- "ExtControl"
      }

    }
  } else {
      pa_message("Using existing Signac/Seurat object")
      if (!inherits(inObj, "Seurat")) {
        stop("Pre-computed input object is not Seurat object!")
      }

      mb <- inObj
      annData <- inObj@meta.data
      if (! (targColumn %in% colnames(annData) ) ) {
        stop("Required annotation column is not available
             in pre-computed input object: ", targColumn)
      }
      if (is.null(ctrlObj)) {
        annInfo <- summary(as.factor(annData[, targColumn]))
        pa_message(paste(capture.output(annInfo)))
        ctrlStatus <- as.numeric(
                  str_split_1(ctrlGrp,pattern = ",") %in% names(annInfo))
        if (any(ctrlStatus == 0) ) {
          stop("Non-tumor control group is not found in annotation: ", ctrlGrp)
        }
      } else {
         pa_message("Merge with external control object")
         mb <- mergeWithControl(mb,ctrlObj,targColumn)
         ctrlGrp <- "ExtControl"
         #print(head(mb@meta.data))
      }
  }

  if (!(dir.exists(resDir))) {
    pa_message(paste0("Creating result directory: ", resDir))
    dir.create(resDir)
  }


  if (!is.null(binSize)) {
    if (is.null(chromLength)) {
      stop("The chromosomes length vector is not provided
              for the bin size adjustment!")
    }
  }

  if (is.null(inObj)) {
    pa_message("Prepare input data...")
    mb <- prepareAnalysis(mb, resDir, sId, annData, targColumn, ctrlObj,
                          performGA, verbose)
  }

  pa_message("Save signal...")
  saveCnvInput(mb, resDir, sId, targColumn)
  if (!is.null(binSize)) {
    pa_message(paste0("Re-format input signal matrix for bin size ", binSize))
    mb <- aggregateBins(mb, resDir, sId, binSize, chromLength)
  }
  #print(head(mb@meta.data))
  if (metaCells) {
    pa_message("Forming meta-cells...")
    pa_message(targColumn)
    if (is.numeric(metaCells)) {
      pa_message("Using custom meta-cell count...")
      metaCount <- metaCells
    } else {
      pa_message("Using default meta-cell count...")
      metaCount <- 5
    }
    extractMetacells(resDir, sId, mb, targColumn,
                     metacell_content = metaCount, verbose )
  }

  pa_message("Write configuration...")
  writeConfig(resDir, sId, ctrlGrp, binSize, metaCells)
  pa_message("Prepared input.")

  invisible(NULL)

}
