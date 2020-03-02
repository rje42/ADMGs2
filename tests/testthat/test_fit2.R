dag20 <- makeGraphChain(20, "directed")

dat <- matrix(combinations(c(2,2)), 4, 20)
dat <- as.data.frame(dat)
names(dat) <- paste("X", 1:20, sep="")

set.seed(143)
dat <- cbind(dat, freq=rpois(4, 50))

test_that("New sparser method works", {
  expect_equal(fitADMG(dat = dat, graph=dag20)$ll, -2807.874671787004)
})

#test_that("New sparser method works for autoFit", {
#  expect_success(autoFit(dat = dat))
#})