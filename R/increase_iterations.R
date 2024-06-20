#' @title Check brmsfit objects save as RDS files
#'
#' @description \code{increase_iterations} checks brmsfit objects save as RDS files.
#' @param path Directory in which to look for .rds files.
#' @param fits Name of the .rds files to be read. Optional.
#' @param check.rhat NOT YET IMPLEMENTED.
#' @param increase NOT YET IMPLEMENTED.
#' @param chains NOT YET IMPLEMENTED.
#' @param cores Number of cores to use for parallelization. Default is 1 (no parallelization).
#' @param pb Logical to control if a progress bar is used. Default is TRUE.
#' @param robust Logical to control the type of central tendency measure as in \code{\link[brms]{summary.brmsfit}}). NOT YET IMPLEMENTED.
#' @param html Logical to control whether results are returned in html format. Useful for creating Rmd or quarto html reports. Is FALSE (default) the table is return as a data frame object. NOT YET IMPLEMENTED.
#' @param verbose Logical to control if messages are printed into the console. NOT YET IMPLEMENTED.
#' @return Returns a data frame with a summary of fitted models. Can be used to make sure all models were run with the same parameters (e.g. before combining models). If \code{html = FALSE} the function will return a data frame, otherwise it will print the estimates in a table in html format. The summary includes: prior, formula, number of iterations, number of chain, thinning, warmup, number of parameters, number of divergent transitions, number of rhats higher than 1.05, tail and bulk effective sample sizes and seed.
#' @export
#' @name increase_iterations
#' @details The function reads all fits saved as rds files in the supplied directory and generates a table listing the parameters used to fit models.
#' @examples
#' {
#' # create directory
#' dir.create(file.path(tempdir(), "rdss"))
#' # run 2 models
#' md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
#' iter = 500, file = file.path(tempdir(), "rdss", "md1"))
#'
#' md2 <- brm(Petal.Length ~ Species, iris, chains = 1,
#' iter = 500, file = file.path(tempdir(), "rdss", "md2"))
#'
#' # check fits
#' increase_iterations(path = file.path(tempdir(), "rdss"))
#' }
#' @seealso \code{\link{extended_summary}}, \code{\link{combine_rds_fits}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: miscellaneous functions to customize brms bayesian regression models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#' }

increase_iterations <- function(path = ".", fits = list.files(path = path, pattern = ".RDS$", ignore.case = TRUE, full.names = TRUE), cores = 1, pb = TRUE, check.rhat = FALSE, increase = "2x", chains = 4, robust = TRUE, html = FALSE, verbose = TRUE){

  # set empty variables to avoid checking errors
  robust <- diverg_transitions <- verbose <- html <- NULL

  if (pb){
    cat("Checking fits (step 1 out of 2):")
  }
  rds_check_df <- check_rds_fits(path = path, cores = cores, fits = fits, pb = pb, html = FALSE)


  # run loop over fits
  fit_table_list <- pblapply_brmsish_int(X = fits, cl = cores, pbar = pb ,function(x){

    # set null objects to avoid errors when checking in CRAN
    Parameter <- NULL

    # read fit
    fit <- try(readRDS(x), silent = TRUE)

    if (!methods::is(fit, 'try-error') & methods::is(fit, 'brmsfit')){

      # get diagnostic quantities
      iterations <- fit$fit@sim$iter
      chains <- posterior::nchains(fit)
      warmup <- fit$fit@sim$warmup
      mod_formula <- as.character(fit$formula[1])
      # diverg_transitions <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
      # percent_transitions <- diverg_transitions / nrow(brms::nuts_params(fit, pars = "divergent__"))
      # diverg_transitions <- paste0(diverg_transitions, " (", percent_transitions, "%)")
      priors <- paste(apply(brms::prior_summary(fit, all = FALSE)[, 2:1], 1, paste, collapse = "-"), collapse = "\n")
      seed <- fit$fit@stan_args[[1]]$seed
      thinning <- fit$fit@stan_args[[1]]$thin
      n_parameters <- length(names(fit$fit@sim$samples[[1]]))

      # get draws
      vars <- posterior::variables(fit)
      fit <- posterior::as_draws_array(fit, variable = vars)

      # get summary
      coef_table <- draw_summary(fit, variables = vars, probs = c(0.025, 0.975), robust = robust)

      # put parameters on a table
      fit_table <-
        data.frame(
          fit = basename(x),
          priors = priors,
          formula = mod_formula,
          iterations = iterations,
          chains = chains,
          thinning = thinning,
          warmup = warmup,
          n_parameters = n_parameters,
          diverg_transitions = diverg_transitions,
          `rhats > 1.05` = sum(coef_table$Rhat > 1.05),
          min_bulk_ESS = min(coef_table$Bulk_ESS),
          min_tail_ESS = min(coef_table$Tail_ESS),
          seed = seed,
          check.names = FALSE
        )

    } else
      fit_table <- NA

    return(fit_table)

  })

  # which were not read
  which_data_frame <- sapply(fit_table_list, is.data.frame)

  not_read <-basename(fits)[!which_data_frame]

  if (length(not_read) > 0 & verbose){

    if (length(not_read) != length(fit_table_list))
      cat("Some RDS files cannot be read (check '.Options$brmsish$not_read_rds')")

    # add not read to options
    options(brmsish = list(not_read_rds = not_read))

    # remove non-data frames
    fit_table_list <- fit_table_list[which_data_frame]
  }

  if (length(fit_table_list) == 0) {
    if(verbose) cat("no readable fits found")
  } else{
    fit_parameters <- do.call(rbind, fit_table_list)

  # remove non brms fits
  fit_parameters <- fit_parameters[!is.na(fit_parameters$formula),]

  if (!html)
    return(fit_parameters) else {

      html_tab <- html_format_fit_table(x = fit_parameters)

      return(html_tab)
    }
  }
}
