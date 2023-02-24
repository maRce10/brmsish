test_that("basic", {

  # run model
  md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
  iter = 1000)

  # print summary
  es <- extended_summary(fit = md1)


  expect_true(is.list(es))
})
