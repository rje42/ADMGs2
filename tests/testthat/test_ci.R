x <- list(list(1,2,3), list(1:3, 4:5, 3:4))

ind1 <- as.ci(x[[1]])
ind2 <- as.ci(x[[2]])

test_that("CIs created correctly", {
  expect_equal(class(ind1), "ci")
  expect_equal(class(ind2), "ci")
  expect_equal(sapply(ind1, class), rep("integer", 3))
  expect_equal(sapply(ind2, class), rep("integer", 3))
})

ind2a <- split_ci(ind2)

test_that("CIs created correctly", {
  expect_equal(sapply(ind2a, class), rep("ci", 6))
  for (i in 1:6) expect_equal(sapply(ind2a[[i]], class), rep("integer", 3))
})

ind2b <- standardize_cis(ind2a)

test_that("CIs created correctly", {
  expect_equal(sapply(ind2b, class), rep("ci", 2))
  for (i in 1:2) expect_equal(sapply(ind2b[[i]], class), rep("integer", 3))
})

