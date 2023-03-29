test_that("findSep works", {
  expect_equal(findSep(gr1, X=3, Y=4), 1:2)
  expect_equal(findSep(gr1, X=1, Y=2), integer(0))
  expect_true(is.na(findSep(gr1, X=1, Y=3)))
})

test_that("listSep works", {
  expect_equal(listSep(gr1, X=3, Y=4), list(1:2))
  expect_equal(listSep(gr1, X=1, Y=2), list(integer(0)))
  expect_equal(listSep(gr1, X=1, Y=3), list())
  
  expect_equal(listSep(gr2, X=1, Y=5), list(2, c(2,4)))
  expect_equal(listSep(gr2, X=1, Y=4), list(2, c(2,5)))
  expect_equal(listSep(gr2, X=1, Y=4:5), list(2))
})

cyc6 <- makeGraphCycle(6)
cyc6b <- makeGraphCycle(6, "bidirected")

test_that("m_connected works", {
  expect_equal(sort.int(m_connected(cyc6, 1, c(3,6))), c(1:3,6))
  expect_equal(sort.int(m_connected(cyc6b, 1, c(3,6))), c(1,2,5,6))
  expect_equal(sort.int(m_connected(cyc6b, 1, 3)), c(1,2,6))
})

