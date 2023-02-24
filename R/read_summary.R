#' @title Read and print the output from 'fit_summary'
#'
#' @description \code{fit_summary} reads and prints the saved output from 'fit_summary'.
#' @usage read_summary(path = ".", fill = "#6DCD59FF", highlight = FALSE)
#' @param path Directory in which to look for the 'fit_summary' output files.
#' @param fill Color for posterior distribution fill. Default is "#6DCD59FF".
#' @param highlight Logical to control if posterior estimates for which the 95\% credible intervals do not overlap with zero are highlighted. Default is FALSE.
#' @return If \code{plot = TRUE} the function reads the saved output of \code{\link{fit_summary}} and prints it.
#' @export
#' @name read_summary
#' @details It prints a summary of brmsfit results that were previously saved by \code{\link{fit_summary}}. Useful when models are too big so reading them many times is to be avoided. It includes a model fit table, a coefficient table and a posterior distribution halfeye graph next to a chain trace plot. Tables are produce in html format so the output is nicely printed when knit in Rmarkdown files. You might have to add 'results = 'as.is', warning = FALSE and adjust fig.width to chunk options in Rmarkdown documents.
#' @examples
#' {
#' # run model
#' md1 <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
#' iter = 500)
#'
#' # save summary
#' fit_summary(fit = md1, save = TRUE, dest.path = tempdir())
#'
#' # read/print summary
#' read_summary(path = file.path(tempdir(), "md1"))
#' }
#' @seealso \code{\link{combine_rds_fits}}, \code{\link{fit_summary}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: random stuff on brms bayesian models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#' }
#'
read_summary <- function(path = ".", fill = "#6DCD59FF", highlight = FALSE){

  # Read fit output
  mod <- readRDS(file.path(path, "fit_table.RDS"))

  # print model formula
  # cat('\n\n')
  cat('\n\n## ', basename(path), '\n\n')

  # print fit summary table
  fit_table <- if (!is.null(mod$fit_table))
    html_format_fit_table(mod$fit_table) else
      html_format_fit_table(mod$model_table)

  print(fit_table)

  # print estimates
  coef_table <- html_format_coef_table(mod$coef_table, fill = fill, highlight = highlight)

  # print fit result table
  print(coef_table)

  # plot
  path <- normalizePath(path)
  cat("![](", file.path(path, "plot.jpeg"), ")", sep = "")
}

