## code to prepare `alstyne` dataset goes here
alstyne <-
  array(c(
    48,  34,  37,  49,  38,  28,  35,  57,
    44,  34,  29,  58,  47,  38,  37,  53,
    1,   5,   7,  11,   3,   8,   5,  18,
    1,   7,   7,   5,   1,   2,   4,  24,
    117, 259, 131, 319, 197, 435, 107, 291,
    111, 253, 131, 320, 202, 392, 103, 294,
    23,  61,  20,  89,  38, 194,  27, 101,
    27,  55,  25,  93,  46, 215,  34, 102
  ), rep(2,6))

alstyne = cbind(expand.grid(list(c(1,0))[rep(1,6)]), c(alstyne))
alstyne <- as.data.frame(alstyne)
names(alstyne) = c("person", "age", "drug", "validation", "success", "prior", "freq")
alstyne$validation <- 1-alstyne$validation


usethis::use_data(alstyne)
