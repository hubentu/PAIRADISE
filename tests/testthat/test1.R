
set.seed(12345)
data("sample_dataset")
pdat <- PDseDataSetFromMat(sample_dataset)
pdat <- pairadise(pdat, numCluster = 2)
res <- results(pdat, sig.level = 0.01, details = TRUE)

test_that("return 2 sig results", {
    expect_equal(nrow(res), 3)
})

test_that("valid counts", {
    expect_equal(ncol(res$latent[1,]$latent), 3)
    expect_equal(ncol(res$latent[3,]$latent), 2)
})
