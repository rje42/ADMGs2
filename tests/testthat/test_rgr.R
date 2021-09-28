set.seed(1903)
admg0 <- rADMG(5, N=10)

dist_n <- lapply(admg0, rADMGdist, dims=2:6, r=FALSE, alpha=2)
dist_n2 <- mapply(probdist, dist_n, graph=admg0)
dist_r <- lapply(admg0, rADMGdist, dims=rep(2,5), r=TRUE, alpha=2)
dist_r2 <- mapply(probdist, dist_r, graph=admg0)
dist_nw <- lapply(admg0, rADMGdist, dims=2:6, r=FALSE, alpha=2, new=TRUE)
dist_nw2 <- mapply(probdist, dist_n, graph=admg0)

test_that("distributions simulated correctly", {
  expect_true(all(dist_n2 >= 0))
  expect_true(all(dist_r2 >= 0))
  expect_true(all(dist_nw2 >= 0))
})
