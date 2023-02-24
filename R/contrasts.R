#' @title Estimate hypothesis testing for all pairwise comparisons
#'
#' @description \code{contrasts} estimates hypothesis testing for all pairwise comparisons of a predictor in brmsfit object.
#' @usage contrasts(fit = NULL, predictor, level.sep = " - ", xlab = "Effect size",
#' gsub.pattern = NULL, gsub.replacement = NULL, n.posterior = 2000,
#' fill = "#6DCD59FF", sort.levels = NULL, html.table  = FALSE,
#' read.file = NULL, plot = FALSE, plot.area.prop = 1, highlight = FALSE)
#' @param fit A brmsfit object.
#' @param predictor The name of the categorical predictor in the model fit for which contrasts will be computed.
#' @param level.sep A character string to separate the levels in the output.
#' @param xlab A character string with the horizontal axis label. Default is "Effect size".
#' @param gsub.pattern A vector with character strings to be replaced
#' @param gsub.replacement A vector with character strings to use for replacement.
#' @param n.posterior Number of posterior samples to use for plotting. Default is 2000.
#' @param fill Color for posterior distribution fill. Default is "#6DCD59FF".
#' @param sort.levels Character vector with the order to be used for levels in the predictor.
#' @param html.table Logical to control whether estimate tables are plotted in html format. Useful for creating Rmd or quarto html reports. Is FALSE (default) the table is return as a data frame object.
#' @param read.file Character string with the name of the .rds file containing the model fit.
#' @param plot Logical to control if posterior distributions of estimates are plotted. Default is FALSE.
#' @param plot.area.prop Positive number to control de proportion of the plotting area of posterior distributions that will be included. Default is 1 (the default area included by \code{\link[ggplot2]{ggplot}}). Useful for adding or removing empty space around distributions.
#' @param highlight Logical to control if posterior estimates for which the 95\% credible intervals do not overlap with zero are highlighted. Default is FALSE.
#' @return If \code{plot = TRUE} the function returns a ggplot object with the posterior distributions of the comparisons between predictor levels. If \code{html = FALSE} the function will return a data frame with estimates for each comparison, otherwise it will print the estimates in a table in html format.
#'
#' @export
#' @name contrasts
#' @details Estimates hypothesis testing for all pairwise comparisons of levels from a categorical predictor. The function \code{\link[brms]{hypothesis}} is used internally.
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
#' @seealso \code{\link{fit_summary}}, \code{\link{read_summary}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: random stuff on brms bayesian models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
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
           highlight = FALSE
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
    levels_df$sign[grep(base_level, levels_df$V1)] <- -1

    names(contrsts) <- paste0(levels_df$V1, level.sep, levels_df$V2)
    contrsts <- gsub(paste0(base_level, " - "), "", contrsts)
    contrsts <- gsub(paste0(" - ", base_level), "", contrsts)
    names(contrsts) <- gsub(predictor, "", names(contrsts))

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

    if (html.table){

      # print estimates
      html_table <- html_format_coef_table(hyp_table, fill = fill,  highlight = highlight)

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

    coef_table2 <- hyp_table
    coef_table2$value <- coef_table2$Estimate
    coef_table2$significance <-
      ifelse(coef_table2$CI_low * coef_table2$CI_high > 0, "sig", "non-sig")
    coef_table2$significance <-
      factor(coef_table2$significance, levels = c("non-sig", "sig"))

    if (highlight)
      col_pointrange <- ifelse(coef_table2$significance == "non-sig", "gray", "black")
    else col_pointrange <- rep("black", nrow(coef_table2))

    fill_values <-
      if (!highlight)
        rep(grDevices::adjustcolor(fill, alpha.f = 0.5), nrow(coef_table2)) else
          if (all(coef_table2$significance == "non-sig"))
            grDevices::adjustcolor(fill, alpha.f = 0.25) else
              if (all(coef_table2$significance == "sig"))
                grDevices::adjustcolor(fill, alpha.f = 0.5) else
                  c(
                    grDevices::adjustcolor(fill, alpha.f = 0.25),
                    grDevices::adjustcolor(fill, alpha.f = 0.5)
                  )

    posteriors$significance <-
      sapply(posteriors$Hypothesis, function(x)
        coef_table2$significance[coef_table2$Hypothesis == x])

    posteriors$sign <-
      sapply(posteriors$Hypothesis, function(x)
        levels_df$sign[levels_df$hypothesis == x])

    posteriors$value <- posteriors$value * posteriors$sign

    posteriors$Hypothesis <- factor(posteriors$Hypothesis, levels = coef_table2$Hypothesis[nrow(coef_table2):1])

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
      ggplot2::geom_point(data = coef_table2) +
      ggplot2::geom_errorbar(data = coef_table2,
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
      ggplot2::labs(x = xlab, y = "Contrasts") +
      ggplot2::xlim(range(c(posteriors$value, 0)) * plot.area.prop)


    print(gg_distribution)
}
    if (!html.table)
    return(hyp_table)

  }
