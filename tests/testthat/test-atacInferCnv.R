test_that("prepare input function parses data correctly", {

  ctrlId = "Normal"
  sId = "MB183_ATAC_test"
  inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz",
                       package = "atacInferCnv")
  sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann.txt",
                     package = "atacInferCnv" )
  resPath = "incorrect_path"
  # no result path
  expect_error(
    prepareAtacInferCnvInput("bad_input", "meta"),
    regexp = "Path to result is not provided!"
  )

  # input mtx incorrect
  expect_error(
    prepareAtacInferCnvInput("bad_input", "meta", resPath),
    regexp = "Input data folder is not found"
  )

  # ann mtx incorrect
  expect_error(
    prepareAtacInferCnvInput(inPath, "meta", resPath),
    regexp = "Annotation file is not found"
  )

  # ann column
  expect_error(
    prepareAtacInferCnvInput(inPath, sAnn, resPath),
    regexp = "Required annotation column is not available"
  )

  # normal cell types
  expect_error(
    prepareAtacInferCnvInput(inPath, sAnn, resPath,
                             targColumn = "cnvBlock", ctrlGrp = "boo"),
    regexp = "Non-tumor control group is not found in annotation"
  )

  expect_error(
    prepareAtacInferCnvInput(resDir= resPath, inObj = "123",
                             targColumn = "cnvBlock"),
    regexp = "Tumor input object is not Seurat"
  )


  expect_error(
    prepareAtacInferCnvInput(inPath, sAnn, resPath,
                             targColumn = "cnvBlock", ctrlObj = 123),
    regexp = "Control input object is not Seurat"
  )


})

test_that("atacInferCnv works correctly on toy data", {

  ctrlId = "Normal"
  sId = "MB183_ATAC_test"
  inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz",
                       package = "atacInferCnv")
  sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann_n30.txt",
                     package = "atacInferCnv" )
  resPath = tempfile()

  expect_message(prepareAtacInferCnvInput(inPath, sAnn, resPath,
                                          targColumn = "cnvBlock",
                                          ctrlGrp = ctrlId,performGA = FALSE),
                 regexp = "Prepared input")

  expect_true(file.exists(file.path(resPath,"sample_raw_counts.txt.gz")))
  expect_true(file.exists(file.path(resPath,"sample_cnv_ref.txt")))
  expect_true(file.exists(file.path(resPath,"sample_cnv_ann.txt")))
  expect_true(file.exists(file.path(resPath,"sample_obj.RDS")))

  # test SingleCellExperiment input
  iObj <- readRDS(file.path(resPath,"sample_obj.RDS"))
  iObj.sce <- as.SingleCellExperiment(iObj)
  resPath2 = tempfile()

  expect_message( prepareAtacInferCnvInput(resDir= resPath2, inObj = iObj.sce,
                      targColumn = "cnvBlock",ctrlGrp = ctrlId, verbose = T),
                  regexp = "Using existing input object")

  expect_true(file.exists(file.path(resPath2,"sample_raw_counts.txt.gz")))
  expect_true(file.exists(file.path(resPath2,"sample_cnv_ref.txt")))
  expect_true(file.exists(file.path(resPath2,"sample_cnv_ann.txt")))


  # test infercnv
  expect_message( runAtacInferCnv(resPath),
                  regexp = "Making the final infercnv heatmap" )
  expect_true(
    file.exists(file.path(resPath,"sample_infercnv","run.final.infercnv_obj")))

  iObj <- readRDS(file.path(resPath,"sample_infercnv","run.final.infercnv_obj"))

  # check params
  expect_true(iObj@options$k_obs_groups == 1)
  expect_true(sum(iObj@options$chr_exclude == c("Y","MT")) == 2)

  # check atac matrix
  expect_true(nrow(iObj@expr.data) == 12952)
  expect_true(ncol(iObj@expr.data) == 30)
  expect_true(max(iObj@expr.data) > 1.5)

  # check pre-defined clustering result
  expect_true(length(iObj@tumor_subclusters$subclusters)  == 3)
  expect_true(length(iObj@tumor_subclusters$subclusters$C1$C1)  == 16)
  expect_true(length(iObj@tumor_subclusters$subclusters$C2$C2)  == 10)
  expect_true(length(iObj@tumor_subclusters$subclusters$Normal$Normal) == 4)

})

