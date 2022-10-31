# takes several brms models stored as RDS files and combine them in a single model. By default the combined model is saved as a RDS file, but it can instead be returned as an object (if \code{save = FALSE}).
# if not supplied the name is obtained from the containing folder
combine_model_rds <-
  function(path,
           dest.path = ".",
           save = TRUE,
           name = NULL,
           overwrite = FALSE,
           suffix = "-combined_model",
           check_data = TRUE) {
    if (is.null(name))
      name <- basename(path)

    file.name <-
      file.path(dest.path, paste0(name, suffix, ".RDS"))

    if (!file.exists(file.name) | overwrite) {
      mods_rds <-
        list.files(
          path = path,
          full.names = TRUE,
          pattern = "\\.RDS$",
          ignore.case = TRUE
        )

      mods_list <- lapply(mods_rds, readRDS)

      mods_comb <- combine_models(mlist = mods_list, check_data = check_data)

      if (save) {
        saveRDS(mods_comb, file.name)
        return(NULL)
      } else
        return(mods_comb)
    }
  }
