data(gr2, package = "MixedGraphs")
int1 <- intrinsicSets(gr2, r=TRUE, sort=3)
int0 <- intrinsicSets(gr2, r=FALSE, sort=3)

test_that("intrinsicSets works OK", {
  expect_equal(int0[order(sapply(int0, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int1[order(sapply(int1, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
})

test_that("intrinsicSets works OK for CADMGs", {
  expect_equal(int0[order(sapply(int0, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int1[order(sapply(int1, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
})

ht1 <- headsTails(gr2, r=TRUE, sort=3)
ht0 <- headsTails(gr2, r=FALSE, sort=3)

test_that("headsTails works OK", {
  expect_equal(ht0$heads, list(1, 2, 3, c(1,3), 2:3, 4, 3:4, 5, c(2,5), c(4,5)))
  expect_equal(ht1$heads, list(1, 2, 3, c(1,3), 2:3, 4, 3:4, 5, c(2,5), c(4,5)))
  expect_equal(ht0$tails, list(2, 4, 5, c(2,4,5), 4:5, integer(0), 5, integer(0), 4, integer(0)))
  expect_equal(ht1$tails, list(2, 4, 5, c(2,4,5), 4:5, integer(0), 5, integer(0), 4, integer(0)))
})
