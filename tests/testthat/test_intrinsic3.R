data(gr2, package = "MixedGraphs")
int1r <- intrinsicSets3(gr2, r=TRUE, sort=3)
int1 <- intrinsicSets3(gr2, r=FALSE, sort=3)
int0r <- intrinsicSets3(grVerma, r=TRUE, sort=3)
int0 <- intrinsicSets3(grVerma, r=FALSE, sort=3)
int5 <- intrinsicSets3(makeGraphChain(5, "b"), r=FALSE, sort=3)
int6 <- intrinsicSets3(makeGraphChain(6, "b"), r=FALSE, sort=3)
int6a <- intrinsicSets3(graphCr("1 -> 2 <-> 3 <-> 4 <-> 5 <- 6"), r=FALSE, sort=3)
int7 <- intrinsicSets3(makeGraphChain(7, "b"), r=FALSE, sort=3)
int7a <- intrinsicSets3(graphCr("1 -> 2 <-> 3 <-> 4 <-> 5 <-> 6 <- 7"), r=FALSE, sort=3)
int6zh <- intrinsicSets(graphCr("1 -> 4 <-> 2 <-> 1 <-> 3 <-> 5 <- 4"), r=FALSE, sort=3)
int6zha <- intrinsicSets3(graphCr("1 -> 4 <-> 2 <-> 1 <-> 3 <-> 5 <- 4"), r=FALSE, sort=3)
int6cyc <- intrinsicSets3(graphCr("2 <-> 1 <-> 3 <-> 4 <-> 6 <-> 5 <-> 2"), r=FALSE, sort=3)


test_that("intrinsicSets3 works OK", {
  expect_equal(int1r, list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int1, list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
  expect_equal(int0r, list(1, 2, 3, 4, c(2,4)))
  expect_equal(int0, list(1, 2, 3, c(2,4)))
  expect_equal(int5, list(1,2,1:2,3,2:3,1:3,4,3:4,2:4,1:4,5,4:5,3:5,2:5,1:5))
  expect_equal(int6, list(1L, 2L, 1:2, 3L, 2:3, 1:3, 4L, 3:4, 2:4, 1:4, 5L, 4:5, 3:5, 
                          2:5, 1:5, 6L, 5:6, 4:6, 3:6, 2:6, 1:6))
  expect_equal(int6a, list(1L, 2L, 3L, 2:3, 4L, 3:4, 2:4, 5L, 4:5, 3:5, 2:5, 6L))
  expect_equal(int7, list(1L, 2L, 1:2, 3L, 2:3, 1:3, 4L, 3:4, 2:4, 1:4, 5L, 4:5, 3:5, 
                          2:5, 1:5, 6L, 5:6, 4:6, 3:6, 2:6, 1:6, 7L, 6:7, 5:7, 4:7, 
                          3:7, 2:7, 1:7))
  expect_equal(int7a, list(1L, 2L, 3L, 2:3, 4L, 3:4, 2:4, 5L, 4:5, 3:5, 2:5, 6L, 5:6, 
                           4:6, 3:6, 2:6, 7L))
  expect_equal(int6zh, list(1L, 2L, 1:2, 3L, c(1L, 3L), 1:3, 4L, c(1L, 2L, 4L), 1:4, 
                            5L, c(1L, 3L, 5L), 1:5))
})

# test_that("intrinsicSets works OK for CADMGs", {
#   expect_equal(int0[order(sapply(int0, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
#   expect_equal(int1[order(sapply(int1, function(x) sum(2^x)))], list(1, 2, 3, 4, 5, 4:5, c(2,4,5), 3:5, 2:5, 1:5))
# })
# 
ht1r <- headsTails3(gr2, r=TRUE, sort=3)
ht1 <- headsTails3(gr2, r=FALSE, sort=3)
ht0r <- headsTails3(grVerma, r=TRUE, sort=3)
ht0 <- headsTails3(grVerma, r=FALSE, sort=3)

test_that("headsTails works OK", {
  expect_equal(ht1r$heads, list(1, 2, 3, c(1,3), 2:3, 4, 3:4, 5, c(2,5), c(4,5)))
  expect_equal(ht1$heads, list(1, 2, 3, c(1,3), 2:3, 4, 3:4, 5, c(2,5), c(4,5)))
  expect_equal(ht1r$tails, list(2, 4, 5, c(2,4,5), 4:5, integer(0), 5, integer(0), 4, integer(0)))
  expect_equal(ht1$tails, list(2, 4, 5, c(2,4,5), 4:5, integer(0), 5, integer(0), 4, integer(0)))
  expect_equal(ht0r$heads, list(1, 2, 3, 4, c(2,4)))
  expect_equal(ht0$heads, list(1, 2, 3, 4))
  expect_equal(ht0r$tails, list(integer(0), 1, 2, 3, c(1,3)))
  expect_equal(ht0$tails, list(integer(0), 1, 2, 1:3))
})
