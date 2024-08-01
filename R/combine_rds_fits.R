#' @title Combine brmsfit objects saved as RDS files.
#'
#' @description \code{combine_rds_fits} checks brmsfit objects saved as RDS files.
#' @usage combine_rds_fits(path, dest.path = ".", save = TRUE, name = NULL, overwrite = FALSE,
#' suffix = "-combined_model", check.data = TRUE, verbose = TRUE, summary = FALSE, ...)
#' @param path Directory in which to look for .rds files. Required.
#' @param dest.path Directory in which to save the combined fit as a RDS file (if \code{save  = TRUE}). Default is the current working directory.
#' @param save Logical to control if the combined fit is saved as a RDS file.
#' @param name Character string with the name of RDS file to be saves (when \code{save  = TRUE})..
#' @param overwrite Logical to control if the RDS file is overwritten (when \code{save  = TRUE}).
#' @param suffix Character string with a suffix to be added to the file name (when \code{save  = TRUE}). Default is "-combined_model".
#' @param check.data Logical. Controls if the data should be checked for being the same across models (as in \code{\link[brms]{brm}}).
#' @param summary Logical to control if a summary of the combined model is also saved. If TRUE \code{\link{extended_summary}} is used.
#' @param verbose Logical to control if messages are printed into the console.
#' @param ... Additional arguments to be passed to \code{\link{extended_summary}} for further customizing summary.
#' @return If \code{save  = TRUE} the combined fit is save as a RDS file, otherwise the function returns a brmsfit object.
#' @export
#' @name combine_rds_fits
#' @details The function takes several brms models stored as RDS files and combine them in a single model. By default the combined model is saved as a RDS file, but it can instead be returned as an object (if \code{save = FALSE}). If not supplied the name is obtained from the containing folder.
#' @examples
#' {
#' # create directory
#' dir.create(file.path(tempdir(), "rdss"))
#' # run 2 models
#' md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
#' iter = 500, file = file.path(tempdir(), "rdss", "md1"))
#'
#' md2 <- brm(Petal.Length ~ Species, iris, chains = 1, iter = 500,
#' file = file.path(tempdir(), "rdss", "md2"))
#'
#' # check fits
#' combine_rds_fits(path = file.path(tempdir(), "rdss"))
#' }
#' @seealso \code{\link{extended_summary}}, \code{\link{check_rds_fits}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: miscellaneous functions to customize brms bayesian regression models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#' }

combine_rds_fits <-
  function(path,
           dest.path = ".",
           save = TRUE,
           name = NULL,
           overwrite = FALSE,
           suffix = "-combined_model",
           check.data = TRUE,
           verbose = TRUE,
           summary = FALSE,
           ...
           ) {

    if (is.null(name))
      name <- basename(path)

    file.name <-
      file.path(dest.path, paste0(name, suffix, ".RDS"))

    if (file.exists(file.name) & !overwrite)
      message(paste0("combined model fit '", basename(file.name), "' already existed (overwrite = FALSE)"))

    if (!file.exists(file.name) | overwrite) {
      mods_rds <-
        list.files(
          path = path,
          full.names = TRUE,
          pattern = "\\.RDS$",
          ignore.case = TRUE
        )

      mods_list <- lapply(mods_rds, readRDS)

    combine_success <- try(mods_comb <- combine_models(mlist = mods_list, check_data = check.data), silent = TRUE)

    if (!is(combine_success, "try-error")){

      if (verbose)
        message(paste0("model fits at '", name, "' successfully combined (", length(mods_rds), " RDS files)"))


              if (save) {
        saveRDS(mods_comb, file.name)

    output <- NULL

      } else {
        output <- mods_comb
        }

      if (summary)
        extended_summary(fit = mods_comb, dest.path = dest.path, save = TRUE, overwrite = overwrite, ...)

    } else
    if (verbose)
      message(paste("model fits at folder", name, "could not be combined"))
    }

  if (!is.null(output))
    return(output)
}
