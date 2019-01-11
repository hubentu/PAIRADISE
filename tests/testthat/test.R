
test_that("Run example", {
    data("sample_dataset")
    pdat <- PDseDataSetFromMat(sample_dataset)
    ## test s4
    expect_s4_class(pdat, "PDseDataSet")
    ## test run
    pdat <- pairadise(pdat, numCluster = 2)
    expect_equal(colnames(rowData(pdat))[3], "outs")
    ## test results
    res <- results(pdat, sig.level = 0.01, details = TRUE)
    expect_equal(nrow(res), 3)
    ## test valid counts
    expect_equal(ncol(res$latent[1,]$latent), 3)
    expect_equal(ncol(res$latent[3,]$latent), 3)
})
