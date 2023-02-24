#' @title Run a brms regression model across a population of phylogenetic trees
#'
#' @description \code{phylogenetic_uncertainty} runs a brms regression model across a population of phylogenetic trees.
#' @usage phylogenetic_uncertainty(formula, data, phylos, sp.id.column, cores = 1,
#' iter = 5000, fit.name, thin = 0, save.fits = FALSE,
#' save.combined = FALSE, path = "./", pb = TRUE, ...)
#' @param formula A model formula. It must not include the random effect term referring to the phylogeny.
#' @param data A data frame containing the data to be used in the model.
#' @param phylos An object of class 'multiPhylo' (see  \code{\link[ape]{multiphylo}}).
#' @param sp.id.column The name of the column containing the 'taxa' label (i.e. the labels that match those in the tree tips).
#' @param cores Number of cores to use for parallelization. Default is 1 (no parallelization).
#' @param iter Integer with the number of the iterations as in \code{\link[brms]{brm}}.
#' @param fit.name Character string with the name of a RDS file to be saved.
#' @param thin Integer with the thinning rate as in \code{\link[brms]{brm}}.
#' @param save.fits Logical to control if single fits are saved as RDS files. Default is FALSE. If 'save.combined' is also FALSE then the model fit is returned into the R environment.
#' @param save.combined Logical to control if the combined single fit is saved as a RDS file. Default is FALSE. If 'save.fits' is also FALSE then the model fit is returned into the R environment.
#' @param path Character string with the directory path in which to save model fit(s) (when either 'save.combined' or 'save.fits' are TRUE). The current working directory is used as default.
#' @param pb Logical to control if a progress bar is used. Default is TRUE.
#' @param ... Additional arguments to be passed to \code{\link[brms]{brm}} for further customizing models.
#' @return The function returns a fit model that combines all submodels with individual phylogenies. Individual submodels can be saved if \code{save.fits = TRUE}.
#' @export
#' @name phylogenetic_uncertainty
#' @details The function allows to take into account phylogenetic uncertainty when running phylogenetically informed regressions by running several models with a population of trees (ideally the highest posterior trees). Individual models are then combined into a single model fit.
#' @examples \dontrun{
#' # example taken from https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
#' phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
#'
#' # quick and dirty trick to make the phylo multiphylo by replicating it
#' phylos <- list(phylo, phylo)
#' class(phylos) <- "multiPhylo"
#'
#' data_simple <- read.table("https://paul-buerkner.github.io/data/data_simple.txt",
#' header = TRUE)
#'
#' # run model
#' pu_mod <- phylogenetic_uncertainty(phen ~ cofactor,                                data = data_simple, sp.id.column = "phylo", phylos = phylos,
#'  iter = 3000, save.fits = FALSE,
#'  save.combined = FALSE, chains = 1,
#'  prior = c(
#'  prior(normal(0, 10), "b"),
#'  prior(normal(0, 50), "Intercept"),
#'  prior(student_t(3, 0, 20), "sd"),
#'  prior(student_t(3, 0, 20), "sigma"))
#'  )
#' }
#' @seealso \code{\link{extended_summary}}, \code{\link{combine_rds_fits}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: miscellaneous functions to customize brms bayesian regression models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#' }


phylogenetic_uncertainty <- function(formula, data, phylos, sp.id.column, cores = 1, iter = 5000, fit.name = "fit", thin = 0, save.fits = FALSE, save.combined = FALSE, path = "./", pb = TRUE, ...){

  if (!save.fits)
    single.fit <- TRUE else single.fit <- FALSE

  if (single.fit & save.combined & file.exists(file.path(path, paste0(fit.name, ".rds"))))
    cat(file.path(path, paste0(fit.name, ".rds")), "file already exists") else {

  vcv.phylos <- lapply(phylos, ape::vcv.phylo)

  # make formula text
  formula <- paste(deparse(formula, width.cutoff = 500), collapse="")

  # add phylo to formula
  formula <- paste0(formula, " + (1|gr(", sp.id.column, ", cov = vcv.phylo))")

  # back to formula format
  formula <- stats::as.formula(formula)

  if (save.fits & !dir.exists(file.path(path, fit.name)))
    dir.create(file.path(path, fit.name))

  # fit preliminary fit
  print("fitting base model (fitted on the first tree)")

  if (save.fits) {
    fits <- list.files(path = file.path(path, fit.name), pattern = fit.name, full.names = TRUE)
    if (length(fits) > 0)
      m.fit <- readRDS(fits[1]) else{
        m.fit <- brm(
          formula = formula,
          data = data,
          data2 = list(vcv.phylo = vcv.phylos[[1]]),
          iter = iter,
          thin = 1,
          file_refit = "always",
          ...
        )

        if (!single.fit)
          saveRDS(object = m.fit, file = file.path(path, fit.name, paste0(fit.name, "-1.rds")))
      }
  } else
    m.fit <- brm(
      formula = formula,
      data = data,
      data2 = list(vcv.phylo = vcv.phylos[[1]]),
      iter = iter,
      thin = 1,
      file_refit = "always",
      ...
    )


  # Loop fit
  print("loop over trees")
  m.fits <- pblapply_brmsish_int(2:length(vcv.phylos), cl = cores, pbar = pb, function(i){

    if (save.fits & !single.fit & !file.exists(file.path(path, fit.name, paste0(fit.name, "-", i, ".rds"))) | !save.fits) {
      up_fit <- stats::update(m.fit,
                       data2 = list(vcv.phylo = vcv.phylos[[i]], iter = iter, thin = thin))

      if (save.fits & !single.fit){
        saveRDS(object = up_fit, file = file.path(path, fit.name, paste0(fit.name, "-", i, ".rds")))

        return(NULL)
      } else return(up_fit)
    }
  })

  # add first fit
  m.fits[[length(m.fits) + 1]] <- m.fit

  if (single.fit){
    no.brmsfit <- which(sapply(m.fits, class) != "brmsfit")

    while(any(sapply(m.fits, class) != "brmsfit")){
      print("fixing failed fits")
      m.fits[no.brmsfit] <- pbapply::pblapply(no.brmsfit, cl = cores, function(i){

        upd.fit <- if (!methods::is(m.fits[[i]], "brmsfit"))
          stats::update(m.fit, data2 = list(vcv.phylo = vcv.phylos[[i]])) else m.fits

        return(upd.fit)
      })
    }

    m.fits <- m.fits[sapply(m.fits, class) == "brmsfit"]

    ### Combine using combine_models()
    m.fits_comb <- combine_models(mlist = m.fits)

    if (save.combined)
      saveRDS(m.fits_comb, file.path(path, paste0(fit.name, ".rds")))
  } else
    m.fits_comb <- NULL
    return(m.fits_comb)
    }
}
