gr_ud1 <- graphCr("x1 - x2", "x1 - x4", "x2 - x5", "x4 - x5", "x1 -> x3", "x2 -> x3", "x4 -> x3", "x5 -> x3")
gr_ud1a <- graphCr("x1 - x2", "x1 - x4", "x2 - x5", "x4 - x5", "x1 -> x3", "x2 -> x3", "x4 -> x3", "x5 -> x3", mode="adjMatrix")
gr_ud1b <- graphCr("x1 - x2", "x1 - x4", "x2 - x5", "x4 - x5", "x1 -> x3", "x2 -> x3", "x4 -> x3", "x5 -> x3", mode="edgeMatrix")
gr_ud1c <- graphCr("x1 - x2", "x1 - x4", "x2 - x5", "x4 - x5", "x1 -> x3", "x2 -> x3", "x4 -> x3", "x5 -> x3", mode="eList")
gr0 <- makeGraphEmpty(3)

ssr <- list(1L, 2L, 1:2, 3L, c(1L, 3L), 2:3, 1:3, 4L, c(1L, 4L), 3:4, 
            c(1L, 3L, 4L), 2:4, 1:4, 5L, c(2L, 5L), c(3L, 5L), c(1L, 3L, 5L), 
            c(2L, 3L, 5L), c(1L, 2L, 3L, 5L), 4:5, 3:5, c(1L, 3L, 4L, 5L), 2:5, 1:5)

test_that("subsetRep works correctly", {
  expect_equal(subsetRep(gr_ud1, sort = 3), ssr)
  expect_equal(subsetRep(gr_ud1a, sort = 3), ssr)
  expect_equal(subsetRep(gr_ud1b, sort = 3), ssr)
  expect_equal(subsetRep(gr_ud1c, sort = 3), ssr)
  expect_equal(subsetRep(gr0, sort = 3), list(1L, 2L, 3L))
})
