test_that("basic", {

  # create directory
  dir.create(file.path(tempdir(), "rdss"))
  # run 2 models
  md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
  iter = 1000, file = file.path(tempdir(), "rdss", "md1"))

  md2 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1, iter = 1000,
  file = file.path(tempdir(), "rdss", "md2"))

  # check fits
  crf <- capture_message(combine_rds_fits(path = file.path(tempdir(), "rdss"), dest.path = file.path(tempdir(), "rdss")))

  unlink(list.files(path = file.path(tempdir(), "rdss"), full.names = TRUE))

  expect_true(crf$message == "model fits at 'rdss' successfully combined (2 RDS files)\n")
})
