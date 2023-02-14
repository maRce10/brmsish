# formula must not include the random effect for species
# the function will read models where it stop last time
# save.models save individual (single tree) models (not the combined one)
# if save.models false then single model = TRUE
phylogenetic_uncertainty <- function(formula, data, phylos, sp.id.column, cores = 1, iter = 5000, model.name, thin = 0, save.models = FALSE, save.combined = FALSE, single.model = FALSE, path = "./", ...){

  if (!save.models)
    single.model <- TRUE

  if (single.model & save.combined & file.exists(file.path(path, paste0(model.name, ".rds"))))
    cat(file.path(path, paste0(model.name, ".rds")), "file already exists") else {

  vcv.phylos <- lapply(phylos, ape::vcv.phylo)

  # make formula text
  formula <- paste(deparse(formula, width.cutoff = 500), collapse="")

  # add phylo to formula
  formula <- paste0(formula, " + (1|gr(", sp.id.column, ", cov = vcv.phylo))")

  # back to formula format
  formula <- as.formula(formula)

  if (save.models & !dir.exists(file.path(path, model.name)))
    dir.create(file.path(path, model.name))

  # fit preliminary model
  print("fit or read base model (fitted on the first tree)")

  if (save.models) {
    models <- list.files(path = file.path(path, model.name), pattern = model.name, full.names = TRUE)
    if (length(models) > 0)
      m.fit <- readRDS(models[1]) else{
        m.fit <- brm(
          formula = formula,
          data = data,
          data2 = list(vcv.phylo = vcv.phylos[[1]]),
          iter = iter,
          thin = 1,
          file_refit = "always",
          ...
        )

        if (!single.model)
          saveRDS(object = m.fit, file = file.path(path, model.name, paste0(model.name, "-1.rds")))
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


  # Loop Model
  print("loop over trees")
  m.fits <- pbapply::pblapply(2:length(vcv.phylos), cl = cores, function(i){

    if (save.models & !single.model & !file.exists(file.path(path, model.name, paste0(model.name, "-", i, ".rds"))) | !save.models) {
      up_fit <- stats::update(m.fit,
                       data2 = list(vcv.phylo = vcv.phylos[[i]], iter = iter, thin = thin))

      if (save.models & !single.model){
        saveRDS(object = up_fit, file = file.path(path, model.name, paste0(model.name, "-", i, ".rds")))

        return(NULL)
      } else return(up_fit)
    }
  })

  # add first model
  m.fits[[length(m.fits) + 1]] <- m.fit

  if (single.model){
    no.brmsfit <- which(sapply(m.fits, class) != "brmsfit")

    while(any(sapply(m.fits, class) != "brmsfit")){
      print("fixing failed models")
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
      saveRDS(m.fits_comb, file.path(path, paste0(model.name, ".rds")))
  } else
    m.fits_comb <- NULL
    return(m.fits_comb)
    }
}
