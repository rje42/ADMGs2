dag1 <- graphCr("1 -> 2 -> 4 <- 3 -> 2")
dag2 <- graphCr("1 -> 2 -> 4 <- 3 <- 1")
dag3 <- graphCr("5 -> 1 -> 2 -> 4 <- 3 <- 1")

p1 <- rDAGdist(10, dag1)
p2 <- rDAGdist(10, dag2)
p3 <- rDAGdist(10, dag3)

test_that("rDAGdist works as expected", {
  expect_true(all(capply(p1, checkCI, 1, 4, 2:3)))
  expect_true(all(capply(p1, checkCI, 1, 3)))
  expect_true(all(capply(p2, checkCI, 1, 4, 2:3)))
  expect_true(all(capply(p2, checkCI, 2, 3, 1)))
  expect_true(all(capply(p3, checkCI, 2:4, 5, 1)))
  expect_true(all(capply(p3, checkCI, 2, 3, c(1,5))))
  expect_true(all(capply(p3, checkCI, 4, 1, 2:3)))
})
