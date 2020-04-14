gr_a <- graphCr("1 -> 2 -> 3 <-> 4 <-- 2, 1 -> 4")
data(grIV, package="MixedGraphs")

ht1 <- headsTails(gr_a, r = FALSE, sort = 3)
ht1a <- headsTails(gr_a, r = FALSE, max_head = 2, sort = 3)

ht2 <- headsTails(grIV, r = FALSE, sort = 3)
ht2a <- headsTails(grIV, r = FALSE, max_head = 1, sort = 3)

ht3 <- headsTails(grVerma, r = FALSE, sort = 3)
ht3a <- headsTails(grVerma, r = FALSE, max_head = 1, sort = 3)
h3b <- unlist(purrr::transpose(headsTails(grVerma, r = FALSE, sort = 3, 
                                          by_district = TRUE))$heads, 
              recursive = FALSE)
h3b <- h3b[order(sapply(h3b, function(x) sum(2^x)))]

test_that("headsTails works for ordinary heads", {
  expect_true(setsetequal(ht1$heads, ht1a$heads))
  expect_true(setsetequal(ht1$tails, ht1a$tails))

  expect_true(setsetequal(ht2$heads, ht2a$heads))
  expect_true(setsetequal(ht2$tails, ht2a$tails))

  expect_true(setsetequal(ht3$heads, ht3a$heads))
  expect_true(setsetequal(ht3$tails, ht3a$tails))
  expect_true(setsetequal(ht3$heads, h3b))
})

### now try recursive heads and tails

ht1r <- headsTails(gr_a, r = TRUE, sort=3)
# ht1a <- headsTails(gr1, r = TRUE, max_head = 2)

ht2r <- headsTails(grIV, r = TRUE, sort=3)
# ht2a <- headsTails(gr2, r = TRUE, max_head = 1)

test_that("headsTails works for recursive heads", {
  expect_true(setsetequal(ht1$heads, ht1r$heads))
  expect_true(setsetequal(ht1$tails, ht1r$tails))

  expect_true(setsetequal(ht2$heads, ht2r$heads))
  expect_true(setsetequal(ht2$tails, ht2r$tails))
})
