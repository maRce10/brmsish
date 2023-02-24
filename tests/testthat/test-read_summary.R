test_that("basic", {

  # run model
  md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
  iter = 1000)

  # save summary
  extended_summary(fit = md1, save = TRUE, dest.path = tempdir())

  # read/print summary
  rs <- read_summary(path = file.path(tempdir(), "md1"))

  expect_true(is.null(rs))
})
