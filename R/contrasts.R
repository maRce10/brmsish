#' @title Estimate hypothesis testing for all pairwise comparisons of categorical predictor levels
#'
#' @description \code{contrasts} estimates hypothesis testing for all pairwise comparisons of the levels of a categorical in a brmsfit object.
#' @usage contrasts(fit = NULL, predictor, level.sep = " - ", xlab = "Effect size",
#' gsub.pattern = NULL, gsub.replacement = NULL, n.posterior = 2000,
#' fill = "#6DCD59FF", sort.levels = NULL, html.table  = FALSE,
#' read.file = NULL, plot = FALSE, plot.area.prop = 1, highlight = FALSE,
#' non.zero = FALSE)
#' @param fit A brmsfit object.
#' @param predictor The name of the categorical predictor in the model fit for which contrasts will be computed. Note that the predictor must have at least 3 levels.
#' @param level.sep A character string to separate the levels in the output.
#' @param xlab A character string with the horizontal axis label. Default is "Effect size".
#' @param gsub.pattern A vector with character strings to be replaced
#' @param gsub.replacement A vector with character strings to use for replacement.
#' @param n.posterior Number of posterior samples to use for plotting. Default is 2000.
#' @param fill Color for posterior distribution fill. Default is "#6DCD59FF".
#' @param sort.levels Character vector with the order to be used for levels in the predictor.
#' @param html.table Logical to control whether estimate tables are plotted in html format. Useful for creating dynamic reports (i.e. Rmd or quarto html reports). Is FALSE (default) the table is return as a data frame object. You might have to add 'results = 'as.is' to chunk options in dynamic reports.
#' @param read.file Character string with the name of the .rds file containing the model fit.
#' @param plot Logical to control if posterior distributions of estimates are plotted. Default is FALSE.
#' @param plot.area.prop Positive number to control de proportion of the plotting area of posterior distributions that will be included. Default is 1 (the default area included by \code{\link[ggplot2]{ggplot}}). Useful for adding or removing empty space around distributions.
#' @param highlight Logical to control if posterior estimates for which the 95\% credible intervals do not overlap with zero are highlighted. Default is FALSE.
#' @param non.zero Logical to determine if predictor level is compared against zero instead. Default is FALSE.
#' @return If \code{plot = TRUE} the function returns a ggplot object with the posterior distributions of the comparisons between predictor levels. If \code{html = FALSE} the function will return a data frame with estimates for each comparison, otherwise it will print the estimates in a table in html format.
#'
#' @export
#' @name contrasts
#' @details Estimates hypothesis testing for all pairwise comparisons of levels from a categorical predictor. The function \code{\link[brms]{hypothesis}} is used internally. Alternatively, if argument \code{non.zero = TRUE} the function evaluates whether each level of the predictor is different from zero.
#'
#' Note that comparisons (i.e. contrasts) of categorical predictor levels when additional predictors are also included in the model are computed at the baseline (categorical predictors) or 0 (continuous predictors) value of the additional predictors. Mean-centering on additional continuous predictors can be used to ensure that the mean value of continuous predictors is used as baseline instead (Schielzeth 2010).
#' @examples
#' {
#' # run model
#' mod <- brm(Petal.Length ~ Species, iris, chains = 1, iter = 500)
#'
#' # compute constrasts with plot
#' contrasts(fit = mod, predictor = "Species", html.table = TRUE, plot = TRUE)
#' # compute constrasts without plot
#' contrasts(fit = mod, predictor = "Species", html.table = TRUE, plot = FALSE)
#' }
#' @seealso \code{\link{extended_summary}}, \code{\link{read_summary}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: miscellaneous functions to customize brms bayesian regression models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#'
#' Schielzeth, H. (2010), Simple means to improve the interpretability of regression coefficients. Methods in Ecology and Evolution, 1: 103-113. https://doi.org/10.1111/j.2041-210X.2010.00012.x
#' }
contrasts <-
  function(fit = NULL,
           predictor,
           level.sep = " - ",
           xlab = "Effect size",
           gsub.pattern = NULL,
           gsub.replacement = NULL,
           n.posterior = 2000,
           fill = "#6DCD59FF",
           sort.levels = NULL,
           html.table  = FALSE,
           read.file = NULL,
           plot = FALSE,
           plot.area.prop = 1,
           highlight = FALSE,
           non.zero = FALSE
           ) {

    if (is.null(fit) & is.null(read.file))
      stop("either 'fit' or 'read.file' must be supplied")

    # object for avoiding errors with ggplot functions when checking package
    significance <-
      value <-
      variable <-
      CI_high <-
      CI_low <-
      Hypothesis <- Parameter <- chain <- iteration <- NULL

    if (is.null(fit) & !is.null(read.file))
      fit <- readRDS(read.file)
    if (!brms::is.brmsfit(fit))
      stop("'fit' is not a 'brmsfit' object")

    variables <- posterior::variables(fit)
    incl_classes <- c(
      "b", "bs", "bcs", "bsp", "bmo", "bme", "bmi", "bm",
      brms:::valid_dpars(fit), "delta", "lncor", "rescor", "ar", "ma", "sderr",
      "cosy", "cortime", "lagsar", "errorsar", "car", "sdcar", "rhocar",
      "sd", "cor", "df", "sds", "sdgp", "lscale", "simo"
    )
    incl_regex <- paste0("^", brms:::regex_or(incl_classes), "(_|$|\\[)")

    variables <- variables[grepl(incl_regex, variables)]
    fit_levels <- gsub("^b_", "", variables)
    fit_levels <- grep(paste0("^Intercept$|", predictor), fit_levels, value = TRUE)

    # fix  baseline level
    original_levels <-
      gsub(" |&", "", paste0(predictor, unique(fit$data[, predictor])))
    base_level <- setdiff(original_levels, fit_levels)

    if(length(original_levels) < 3)
      stop2("The predictor must have at least 3 levels")

    # get levels
    pred_levels <- as.character(unique(fit$data[, predictor]))

    # sort
    if (!is.null(sort.levels))
      pred_levels <- pred_levels[match(sort.levels, pred_levels)]

    # add predictor name
    pred_levels <- paste0(predictor, pred_levels)

    # create data frame with level pairs
    levels_df <- as.data.frame(t(utils::combn(pred_levels, 2)))

    # remove spaces and &
    levels_df$V1.nospace <- gsub(" |&", "", levels_df$V1)
    levels_df$V2.nospace <- gsub(" |&", "", levels_df$V2)

    # create contrasts in brms syntax
    contrsts <- paste(apply(levels_df[, c("V1.nospace", "V2.nospace")], 1, paste, collapse = " - "), "= 0")

   # convert magnitude for those compare against baseline
    levels_df$sign <- 1

    if (!non.zero)
      levels_df$sign[grep(base_level, levels_df$V1)] <- -1

    names(contrsts) <- paste0(levels_df$V1, level.sep, levels_df$V2)
    contrsts <- gsub(paste0(base_level, " - "), "", contrsts)
    contrsts <- gsub(paste0(" - ", base_level), "", contrsts)
    names(contrsts) <- gsub(predictor, "", names(contrsts))

    if (non.zero) {# compare all against 0
      contrsts <- paste(pred_levels , "+ Intercept = 0")

      names(contrsts) <-  paste(pred_levels , "= 0")

      contrsts[grep(base_level, contrsts)] <- "Intercept = 0"
    }

    levels_df$hypothesis <- names(contrsts)

    if (!is.null(gsub.pattern) & !is.null(gsub.replacement)) {
      if (length(gsub.pattern) != length(gsub.replacement))
        stop2("'gsub.replacement' and 'gsub.pattern' must have the same length")

      for (i in 1:length(gsub.pattern))
        names(contrsts) <-
          gsub(gsub.pattern[i], gsub.replacement[i], names(contrsts))
    }

   # evaluate hypothesis
    hyps <- brms::hypothesis(fit, contrsts)

    hyp_table <-
      hyps$hypothesis[, c("Hypothesis",
                          "Estimate",
                          "Est.Error",
                          "CI.Lower",
                          "CI.Upper")]


    hyp_table$Estimate <- hyp_table$Estimate * levels_df$sign
    hyp_table$`l-95% CI` <- hyp_table$CI.Lower * levels_df$sign
    hyp_table$`u-95% CI` <- hyp_table$CI.Upper * levels_df$sign
    hyp_table$CI.Lower <- hyp_table$CI.Upper <- NULL

    # duplicate table to modify columns for output
    temp_table <- hyp_table

    if (non.zero)
      temp_table$Hypothesis <- gsub(predictor, "", temp_table$Hypothesis) else {

        temp_table$Contrasts <- temp_table$Hypothesis
        temp_table$Hypothesis <- NULL
        temp_table <- temp_table[, c("Contrasts", "Estimate", "Est.Error", "l-95% CI", "u-95% CI")]
      }

    if (html.table){

      # print estimates
      html_table <- html_format_coef_table(temp_table, fill = fill,  highlight = highlight)

      # print fit result table
      print(html_table)
    }

    if (plot){
    # subsample posteriors
    xdrws <- brms::as_draws(hyps$samples)

    names(xdrws)[1:length(contrsts)] <- names(contrsts)

    # only apply thinning if length of posterior < n.posterior
    xdrws <-
      posterior::subset_draws(x = xdrws, variable = names(contrsts))

    if (round(length(xdrws[[1]][[1]]) / n.posterior, 0) >= 2)
      xdrws <-
      posterior::thin_draws(xdrws, round(length(xdrws[[1]][[1]]) / n.posterior, 0))

    merged_xdrws <- posterior::merge_chains(xdrws)
    sub_posts <- as.data.frame(merged_xdrws)[, names(contrsts)]
    names(sub_posts) <- names(contrsts)

    hyp_table$CI_low <- round(hyp_table$`l-95% CI`, digits = 3)
    hyp_table$CI_high <- round(hyp_table$`u-95% CI`, digits = 3)
    hyp_table$`l-95% CI` <- hyp_table$`u-95% CI` <- NULL

    out <-
      lapply(names(contrsts), function(y)
        data.frame(
          variable = y,
          value = sort(sub_posts[, colnames(sub_posts) == y], decreasing = FALSE)
        ))

    posteriors <- do.call(rbind, out)
    posteriors$Hypothesis <-
      factor(posteriors$variable, levels = sort(unique(posteriors$variable)))

    hyp_table$value <- hyp_table$Estimate
    hyp_table$significance <-
      ifelse(hyp_table$CI_low * hyp_table$CI_high > 0, "sig", "non-sig")
    hyp_table$significance <-
      factor(hyp_table$significance, levels = c("non-sig", "sig"))

    if (highlight)
      col_pointrange <- ifelse(hyp_table$significance == "non-sig", "gray", "black")  else col_pointrange <- rep("black", nrow(hyp_table))

    fill_values <-
      if (!highlight)
        rep(grDevices::adjustcolor(fill, alpha.f = 0.5), nrow(hyp_table)) else
          if (all(hyp_table$significance == "non-sig"))
            grDevices::adjustcolor(fill, alpha.f = 0.25) else
              if (all(hyp_table$significance == "sig"))
                grDevices::adjustcolor(fill, alpha.f = 0.5) else
                  c(
                    grDevices::adjustcolor(fill, alpha.f = 0.25),
                    grDevices::adjustcolor(fill, alpha.f = 0.5)
                  )

    posteriors$significance <-
      sapply(posteriors$Hypothesis, function(x)
        hyp_table$significance[hyp_table$Hypothesis == x])

    posteriors$sign <- if (!non.zero)
      sapply(posteriors$Hypothesis, function(x)
        levels_df$sign[levels_df$hypothesis == x]) else 1

    posteriors$value <- posteriors$value * posteriors$sign

    posteriors$Hypothesis <- factor(posteriors$Hypothesis, levels = hyp_table$Hypothesis[nrow(hyp_table):1])

    # fix level labels for plot
    if (non.zero) {
      levels(posteriors$Hypothesis) <- gsub(predictor, "", as.character(levels(posteriors$Hypothesis)))
      hyp_table$Hypothesis <- gsub(predictor, "", as.character(hyp_table$Hypothesis))
    }

    gg_distribution <-
      ggplot2::ggplot(data = posteriors,
                      ggplot2::aes(y = Hypothesis, x = value, fill = significance)) +
      ggplot2::geom_vline(xintercept = 0,
                          col = "black",
                          lty = 2) +
      ggdist::stat_halfeye(
        ggplot2::aes(x = value),
        .width = c(.95),
        normalize = "panels",
        color = "transparent"
      ) +
      ggplot2::scale_fill_manual(values = fill_values, guide = 'none')  +
      ggplot2::geom_point(data = hyp_table) +
      ggplot2::geom_errorbar(data = hyp_table,
                             ggplot2::aes(xmin = CI_low, xmax = CI_high),
                             width = 0) +
      # ggplot2::scale_color_manual(values = col_pointrange) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.ticks.length = ggplot2::unit(0, "pt"),
        plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = xlab, y = if(non.zero) "Hypothesis" else "Contrasts") +
      ggplot2::xlim(range(c(posteriors$value, 0)) * plot.area.prop)


    print(gg_distribution)
}
    if (!html.table)
    return(temp_table)

  }
