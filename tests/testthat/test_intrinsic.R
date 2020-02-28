int1 <- intrinsicSets(gr2, r=TRUE, sort=2)
int0 <- intrinsicSets(gr2, r=FALSE, sort=2)

test_that("intrinsicSets works OK", {
  expect_equal(int0[order(sapply(int0, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int1[order(sapply(int1, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
})

test_that("intrinsicSets works OK for CADMGs", {
  expect_equal(int0[order(sapply(int0, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int1[order(sapply(int1, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
})