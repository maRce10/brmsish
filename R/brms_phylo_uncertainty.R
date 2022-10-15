brms_phylo_uncertainty <- function(phylos, data, formula, parallel, iter = 5000, model.name, thin = 0, save.model = FALSE, single.model = FALSE, path = "./", ...){

  As <- lapply(phylos, ape::vcv.phylo)

  if (save.model & !dir.exists(file.path(path, model.name)))
    dir.create(file.path(path, model.name))

  # fit preliminary model
  print("fit or read base model (fitted on the first tree)")

  if (save.model) {
    models <- list.files(path = file.path(path, model.name), pattern = model.name, full.names = TRUE)
    if (length(models) > 0)
      m.fit <- readRDS(models[1]) else{
        m.fit <- brm(
          formula = formula,
          data = data,
          data2 = list(A = As[[1]]),
          iter = iter,
          thin = 1,
          file_refit = "always",
          ...
        )

        if (!single.model)
          saveRDS(object = m.fit, file = file.path(path, model.name, paste0(model.name, "-1.RDS")))
      }
  } else
    m.fit <- brm(
      formula = formula,
      data = data,
      data2 = list(A = As[[1]]),
      iter = iter,
      thin = 1,
      file_refit = "always",
      ...
    )


  # Loop Model
  print("loop over trees")
  m.fits <- pbapply::pblapply(2:length(As), cl = parallel, function(i){

    if (save.model & !single.model & !file.exists(file.path(path, model.name, paste0(model.name, "-", i, ".RDS"))) | !save.model) {
      up_fit <- stats::update(m.fit,
                       data2 = list(A = As[[i]], iter = iter, thin = thin))

      if (save.model & !single.model){
        saveRDS(object = up_fit, file = file.path(path, model.name, paste0(model.name, "-", i, ".RDS")))

        return(NULL)
      } else return(up_fit)
    }
  })

  if (single.model){
    no.brmsfit <- which(sapply(m.fits, class) != "brmsfit")


    while(any(sapply(m.fits, class) != "brmsfit")){
      print("fixing failed models")
      m.fits[no.brmsfit] <- pbapply::pblapply(no.brmsfit, cl = parallel, function(i){

        upd.fit <- if (!methods::is(m.fits[[i]], "brmsfit"))
          stats::update(m.fit, data2 = list(A = As[[i]])) else m.fits

        return(upd.fit)
      })


    }

    m.fits <- m.fits[sapply(m.fits, class) == "brmsfit"]

    ### Combine using combine_models()
    m.fits_comb <- combine_models(mlist = m.fits)


    if (save.model)
      saveRDS(m.fits_comb, paste0(model.name, ".RDS")) else
        return(m.fits_comb)
  } else return(NULL)

}
