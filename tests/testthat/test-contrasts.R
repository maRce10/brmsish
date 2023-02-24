test_that("basic", {
  # run model
  mod <- brm(Petal.Length ~ Species, iris, chains = 1, iter = 1000)

  # compute constrasts without plot
  a <- contrasts(fit = mod, predictor = "Species", html.table = FALSE, plot = TRUE)

  expect_true(nrow(a) == 3)
})
