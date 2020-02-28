gr0 <- graphCr("1 <-> 2 <-> 3")

test_that("test head orders", {
  expect_equal(headOrder(c(1,2), c(2,3), gr0), 0L)
  expect_equal(headOrder(c(1,2,3), c(2,3), gr0), -1L)
  expect_equal(headOrder(3, 1:3, gr0), 1L)
})