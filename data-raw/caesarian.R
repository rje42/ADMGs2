## code to prepare `caesarian` dataset goes here
caesarian = structure(list(cesarean = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,
                                        1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,
                                        1L, 0L, 1L, 0L, 1L, 0L, 1L), 
                           monitor = c(0L, 0L, 1L, 1L, 0L,
                                       0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L,
                                       0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L), 
                           arrest = c(0L, 0L,
                                      0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L,
                                      0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L), 
                           breech = c(0L,
                                      0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L,
                                      0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
                           nullipar = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                        0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L), 
                           freq = c(3886L, 32L, 2747L, 67L,
                                    116L, 13L, 194L, 77L, 114L, 34L, 88L, 17L, 7L, 4L, 5L, 13L,
                                    2390L, 140L, 2822L, 187L, 173L, 89L, 592L, 371L, 86L, 77L,
                                    61L, 28L, 4L, 19L, 11L, 20L)), 
                      .Names = c("cesarean", "monitor",
                                 "arrest", "breech", "nullipar", "freq"), 
                      row.names = c(25L, 9L,
                                    26L, 10L, 27L, 11L, 28L, 12L, 29L, 13L, 30L, 14L, 31L, 15L, 32L,
                                    16L, 17L, 1L, 18L, 2L, 19L, 3L, 20L, 4L, 21L, 5L, 22L, 6L, 23L,
                                    7L, 24L, 8L), class = "data.frame")


usethis::use_data(caesarian)
