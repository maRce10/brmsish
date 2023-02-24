test_that("basic", {

  phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")

  # quick and dirty trick to make the phylo multphylo by replicating it
  phylos <- list(phylo, phylo)
  class(phylos) <- "multiPhylo"

  data_simple <- read.table(
    "https://paul-buerkner.github.io/data/data_simple.txt",
    header = TRUE
  )

  pu_mod <- phylogenetic_uncertainty(phen ~ cofactor,  pb = FALSE,                                         data = data_simple, sp.id.column = "phylo", phylos = phylos,
                                        iter = 3000, save.fits = FALSE,
                                         save.combined = FALSE, chains = 1,
                                     prior = c(
                                       prior(normal(0, 2), "b"),
                                       prior(normal(0, 10), "Intercept"),
                                       prior(student_t(3, 0, 20), "sd"),
                                       prior(student_t(3, 0, 20), "sigma")
                                      )
                                     )

  expect_true(is.brmsfit(pu_mod))
})
