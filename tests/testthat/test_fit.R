admg0 <- mixedgraph(4)
admg1 <- graphCr("1 -> 2 -> 4 <-> 3 <-> 2")
admg2 <- makeGraphCycle(4, "bidirected")
admg3 <- makeGraphComplete(4, "bidirected")

set.seed(1902)
dat <- rpois(16, 50)
dim(dat) <- rep(2,4)

test_that("fit is as expected", {
  expect_equal(fitADMG(dat, admg0)$ll, -2142.163400)
  expect_equal(fitADMG(dat, admg1)$ll, -2136.827979)
})

fit2 <- fitADMG(dat, admg2)

test_that("fit is as expected 2", {
  expect_equal(fit2$ll, -2136.91013111)
  expect_equal(class(fit2), "mixed_fit")
  expect_equal(class(summary(fit2)), "mixed_fit_summary")
})

data(twins)
fit30 <- fitADMG(twins, admg3)
fit3 <- fitADMG(twins, admg2)

test_that("fit is as expected 3", {
  expect_equal(2*(fit30$ll - fit3$ll), 28.8833252)
})

set.seed(1902)
dat2 <- rpois(81, 50)
dim(dat2) <- rep(3,4)
fit4 <- fitADMG(dat2, admg3, tol=1e-10)

test_that("fit is as expected 4", {
  expect_equal(fit4$ll, sum(dat2*log(dat2/sum(dat2))))
})
