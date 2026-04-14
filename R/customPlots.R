
#' Function to plot CNV blocks
#'
#' This function creates a plot for CNV assigned/identified subclones
#' @param resDir Path to the result directory with input
#' @param infercnvObj InferCnv result object. Default: NULL
#' If NULL, then the object will be loaded from result directory.
#' @param save Set TRUE to save plot to the result directory. Default: FALSE
#' @param verbose Detailed output, progress messages, default TRUE
#' @return Returns a list of plots
#' @examples
#' resPath = tempfile()
#' inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz",
#'                       package = "atacInferCnv")
#' sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann_n30.txt",
#'                       package = "atacInferCnv")
#' prepareAtacInferCnvInput(inPath,sAnn,resPath, targColumn = "cnvBlock",
#'                         ctrlGrp = "Normal", performGA = FALSE)
#' runAtacInferCnv(resPath)
#' plotCnvBlocks(resPath)
#' @export

plotCnvBlocks <- function( resDir, infercnvObj = NULL,
                            save = FALSE, verbose = TRUE) {

  cnvDir <- file.path(resDir,"sample_infercnv")
  if (!(dir.exists(cnvDir))) {
    stop("The result directory with InferCNV input does not exist:",resDir)
  }

  if (verbose) {
    message("Load InferCNV result...")
  }

  if (is.null(infercnvObj)) {
    infercnv_obj <- readRDS(file.path(cnvDir, "run.final.infercnv_obj"))
  } else {
    if (!inherits(infercnvObj, "infercnv")) {
      stop("The input object is not infercnv  object!")
    }
    infercnv_obj <- infercnvObj
  }

  obsFile <- file.path(cnvDir,"infercnv.observations.txt")
  if (file.exists(obsFile)) {
    cnvMtx <- read.table(obsFile,check.names = FALSE)
  } else {
    warning("Output file infercnv.observations.txt does not exist in result
           folder, using expr matrix.\n
         Update InferCnv version (>= 1.3.3 ) to support normalized output.")
    cIds <- unlist(infercnv_obj@observation_grouped_cell_indices)
    cnvMtx <- infercnv_obj@expr.data[ ,cIds]
  }
  gnMtx <- infercnv_obj@gene_order
  chrVec <- as.character(gnMtx[, 1])
  cMax <- max(rowMeans(cnvMtx)) * 1.2 #max(cnvMtx)
  cMin <- min(rowMeans(cnvMtx)) * 0.8 #min(cnvMtx)
  cMean <- mean(rowMeans(cnvMtx))

  # draw lines
  curChr <- "1"
  breaks <- c(1)
  chrNames <- c(curChr)
  for (i in seq_len(nrow(gnMtx))) {
    chrId <- as.character(gnMtx[i,1] )
    if (chrId != curChr) {
      #print(chrId)
      breaks <- c(breaks,i)
      curChr <- chrId
      chrNames <- c(chrNames,chrId)
    }
  }


  if (infercnv_obj@options$k_obs_groups  == 1) {
    blocks <- c("Full", names(infercnv_obj@tumor_subclusters$subclusters))
  } else {
    hcBlocks <- cutree(infercnv_obj@tumor_subclusters$hc$all_observations,
                       infercnv_obj@options$k_obs_groups )
    hcBlocksAdj <- paste0("C",hcBlocks)
    names(hcBlocksAdj) <- names(hcBlocks)
    blocks <- c("Full", unique(hcBlocksAdj))
  }

  resPlots <- list()
  for (targ in blocks) {
    if (verbose) {
      message(targ)
    }
    if (targ == "Full") {
      vals<-rowMeans(cnvMtx)
    } else {
      if (infercnv_obj@options$k_obs_groups  == 1) {
        cellIds<-names(infercnv_obj@tumor_subclusters$subclusters[[targ]][[1]])
        if (sum(cellIds %in% colnames(cnvMtx)) == 0) {
          if (verbose) {
            message("Skipping...")
          }
          next
        }
      } else {
        cellIds <- names(hcBlocksAdj)[hcBlocksAdj == targ]
      }
      vals <- rowMeans(cnvMtx[,cellIds])
    }

    plotDf <- data.frame(
      x = seq_along(vals), value = vals, chr = chrVec,
      stringsAsFactors = FALSE
    )

    p <- ggplot(plotDf, aes(x = x, y = value, color = value)) +
      geom_point(size = 0.2) +
      scale_color_gradient(
        low = "red", high = "green",  guide = "none" #, midpoint = 1, # scg2
      ) +
    geom_hline(yintercept = 1, color = "grey50", linetype = "dashed") +
    labs(
      title = paste0("inferCNV signal: ", targ),
      x = "Genomic position",
      y = "Modified signal"
    ) +
    coord_cartesian(ylim = c(cMin, cMax)) +
    scale_x_continuous(breaks=breaks,labels =  chrNames) +
    theme_bw() +
    theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
     # axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

    if (length(breaks) > 0) {
      for (i in seq_len(breaks)) {
        p <- p +
          geom_vline(xintercept = breaks[i],
                     linetype = "dashed", color = "black")

      }
    }

    resPlots[[targ]] <- p

    if (save) {
      resName <- file.path(cnvDir,paste0("subclone_",targ,"_CNV_plot.pdf"))
      ggsave(filename = resName, plot = p, width = 14, height = 6)
    }

  }


   return(resPlots)
}
