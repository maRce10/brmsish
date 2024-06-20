#' @title Print a summary of brmsfit results
#'
#' @description \code{extended_summary} prints a summary of brmsfit results.
#' @param fit A brmsfit object.
#' @param draws A draws object.
#' @param gsub.pattern A vector with character strings to be replaced
#' @param gsub.replacement A vector with character strings to use for replacement.
#' @param xlab A character string with the horizontal axis label. Default is "Effect size".
#' @param ylab A character string with the vertical axis label. Default is "Parameter".
#' @param n.posterior Number of posterior samples to use for plotting. Default is 2000.
#' @param fit.name Character string to be added as title. If not supplied and 'read.file' is supplied the name is obtained from the file.
#' @param read.file Character string with the name of a RDS file to be read, if supplied 'fit' is ignored.
#' @param read.file Character string with the name of the .rds file containing the model fit.
#' @param plot.area.prop Positive number to control de proportion of the plotting area of posterior distributions that will be included. Default is 1 (the default area included by \code{\link[ggplot2]{ggplot}}). Useful for adding or removing empty space around distributions.
#' @param remove.intercepts Logical to control if intercepts are included in the output.
#' @param fill Color for posterior distribution fill. Default is "#6DCD59FF".
#' @param trace.palette Color palette function for trace plotting. Default is \code{\link[viridis]{viridis}}.
#' @param effects Character vector with the name of the effects (e.g. predictor variables) to be included. Optional. Note that if effect names have been modified using 'gsub.pattern' and 'gsub.replacement' the modified names should be used instead. Effects in plots and tables will also be sorted as in 'effects' so this can be used to customize order.
#' @param save Logical to control if the summary is saved instead of printed. If TRUE a RDS file with the function output (a list) and a jpeg file is saved into 'dest.path'.
#' @param dest.path Directory in which to save results (if \code{save = TRUE}). The current working directory is used as default.
#' @param overwrite Logical to control if saved results are overwritten. Defaul is FALSE.
#' @param robust Logical to control the type of central tendency measure as in \code{\link[brms]{summary.brmsfit}}).
#' @param width Width of posterior distribution + trace plots as in \code{\link[ggplot2]{ggsave}}). Default is 8 (in).
#' @param height Height of posterior distribution + trace plots as in \code{\link[ggplot2]{ggsave}}). Default is 'dynamic' which means that the height will increase as more panels (predictors) are added.
#' @param highlight Logical to control if posterior estimates for which the 95\% credible intervals do not overlap with zero are highlighted. Default is FALSE.
#' @param print.name Logical to control if the name of the model fit is printed (when \code{plot = TRUE}).
#' @param trace Logical to control if chain trace plots are generated. Default is \code{TRUE}.
#' @param return Logical to control if tables and plots are returned as a list instead of printed. Default is \code{FALSE} (printing tables and plots).
#' @param spread.type Character string with the name of the method to quantify the spread of posterior distributions. Two options are available: "HPDI" (highest posterior density interval) and "MAD" (mean absolute deviation, default).
#' @return If \code{plot = TRUE} the function returns a model fit table, a coefficient table and a posterior distribution halfeye graph. If \code{save = TRUE} this objects are saved as a RDS file along with jpeg file.
#' @export
#' @name extended_summary
#' @details It prints a summary of brmsfit results. It includes a model fit table, a coefficient table and a posterior distribution halfeye graph next to a chain trace plot. Tables are produce in html format so the output is nicely printed when knit in Rmarkdown files. You might have to add 'results = 'as.is' to chunk options in Rmarkdown documents.
#' @examples
#' {
#' # run model
#' mod <- brm(Petal.Length ~ Petal.Width + Species, iris, chains = 1,
#' iter = 500, file = file.path(tempdir(), "mod"))
#'
#' # print from fit
#' extended_summary(fit= mod)
#'
#' # print summary reading the rds file
#' extended_summary(read.file = file.path(tempdir(), "mod.rds"))
#' }
#' @seealso \code{\link{combine_rds_fits}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas (2022), brmsish: miscellaneous functions to customize brms bayesian regression models. R package version 1.0.0.
#'
#' Paul-Christian Buerkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
#' }

extended_summary <-
  function(fit = NULL,
           draws = NULL,
           gsub.pattern = NULL,
           gsub.replacement = NULL,
           xlab = "Effect size",
           ylab = "Parameter",
           n.posterior = 2000,
           fit.name = NULL,
           read.file = NULL,
           plot.area.prop = 1,
           remove.intercepts = FALSE,
           fill = "#6DCD59FF",
           trace.palette = viridis::viridis,
           effects = NULL,
           save = FALSE,
           dest.path = ".",
           overwrite = FALSE,
           robust = FALSE,
           width = 8,
           height = "dynamic",
           highlight = FALSE,
           print.name = TRUE,
           trace = TRUE,
           return = FALSE,
           spread.type = "MAD"
           ) {

     # create objects just to avoid errors with ggplot functions when checking package
      `l-95% CI` <-
      `u-95% CI` <-
      significance <-
      value <-
      variable <-
      CI_high <-
      CI_low <-
      Hypothesis <- Parameter <- chain <- iteration <- pa_comb_mod <- NULL

    if (is.null(fit) & is.null(read.file) & is.null(draws))
      stop("either 'fit', 'draws' or 'read.file' must be supplied")

    if (!is.null(fit) & is.null(fit.name))
      fit.name <- deparse(substitute(fit)) else
        if (!is.null(read.file) & is.null(fit.name))
          fit.name <- gsub("\\.rds$", "", basename(read.file), ignore.case = TRUE)

    # skip everything if save TRUE and folder exists
    if (!dir.exists(file.path(dest.path, fit.name)) & save | !save | save & overwrite){

    if (is.null(fit) & !is.null(read.file))
      fit <- readRDS(read.file)

    if (!is.null(fit))
      variables <- posterior::variables(fit)
    # if (!is.null(draws))
    #   variables <- grep("b_", names(draws), fixed = TRUE, value = TRUE)


    incl_classes <- c(
      "b", "bs", "bcs", "bsp", "bmo", "bme", "bmi", "bm",
      brms:::valid_dpars(fit), "delta", "lncor", "rescor", "ar", "ma", "sderr",
      "cosy", "cortime", "lagsar", "errorsar", "car", "sdcar", "rhocar",
      "sd", "cor", "df", "sds", "sdgp", "lscale", "simo"
    )
    incl_regex <- paste0("^", brms:::regex_or(incl_classes), "(_|$|\\[)")
    variables <- variables[grepl(incl_regex, variables)]
    betas <- grep("^b_", variables, value = TRUE)

    # remove intercept betas
    if (remove.intercepts)
      betas <- grep("b_Intercept", betas, value = TRUE, invert = TRUE)

    iterations <- fit$fit@sim$iter
    chains <- posterior::nchains(fit)
    warmup <- fit$fit@sim$warmup
    mod_formula <- as.character(fit$formula[1])
    diverg_transitions <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
    percent_transitions <- diverg_transitions / nrow(brms::nuts_params(fit, pars = "divergent__"))
    diverg_transitions <- paste0(diverg_transitions, " (", percent_transitions, "%)")
    priors <- paste(apply(brms::prior_summary(fit, all = FALSE)[, 2:1], 1, paste, collapse = "-"), collapse = "\n")
    seed <- fit$fit@stan_args[[1]]$seed
    thinning <- fit$fit@stan_args[[1]]$thin


    # replace fit with draws (to avoid having several huge objects)
    fit <- posterior::as_draws_array(fit, variable = betas)

    coef_table <- draw_summary(draws = fit, variables = betas, probs = c(0.025, 0.975), robust = robust, spread.type = spread.type)

    min_bulk_ESS <- min(coef_table$Bulk_ESS)
    min_tail_ESS <- min(coef_table$Tail_ESS)
    rhats_1.05 <- sum(coef_table$Rhat > 1.05)

    if (!is.null(draws)){

      tidy_coef_table <-
        tidybayes::summarise_draws(pa_comb_mod,
                                   median,
                                   ~ quantile(.x, probs = c(0.025, 0.975)),
                                   default_convergence_measures())

      min_bulk_ESS <- min(tidy_coef_table$ess_bulk)
      min_tail_ESS <- min(tidy_coef_table$ess_tail)
      rhats_1.05 <- sum(tidy_coef_table$rhat > 1.05)

      }

    fit_table <-
      data.frame(
        priors = priors,
        formula = mod_formula,
        iterations = iterations,
        chains = chains,
        thinning = thinning,
        warmup = warmup,
        diverg_transitions = diverg_transitions,
        `rhats > 1.05` = rhats_1.05,
        min_bulk_ESS = min_bulk_ESS,
        min_tail_ESS = min_tail_ESS,
        seed = seed
      )

    # thin before getting chain data for trace plot
    if (round(fit_table$iterations / n.posterior, 0) >= 2)
      fit <-
      posterior::thin_draws(fit, fit_table$iterations / n.posterior, 0)

    # data by chain for trace plot
    posteriors_by_chain <- do.call(rbind, lapply(betas, function(x) {

      X <- as.data.frame(fit[,, x])
      names(X) <- paste("chain", 1:ncol(X))
      X <- stack(X)
      X$variable <- x
      X$iteration <-
        round(seq(1, fit_table$iterations, length.out = nrow(X) / fit_table$chains))
      names(X) <- c("value", "chain", "variable", "iteration")
      return(X)
    }))

    # merge chains
    fit <- as.data.frame(posterior::merge_chains(fit))
    names(fit) <- betas

    fit <- do.call(rbind, lapply(betas, function(y)
      data.frame(
        variable = y,
        value = sort(fit[, colnames(fit) == y], decreasing = FALSE)
      )))

    fit$variable <-
      factor(fit$variable, levels = sort(unique(fit$variable)))

    coef_table2 <- coef_table
    coef_table2$variable <-
      factor(rownames(coef_table2))
    coef_table2$value <- coef_table2$Estimate
    coef_table2$significance <-
      ifelse(coef_table2$`l-95% CI` * coef_table2$`u-95% CI` > 0, "sig", "non-sig")
    coef_table2$significance <-
      factor(coef_table2$significance, levels = c("non-sig", "sig"))


    if (!is.null(gsub.pattern) & !is.null(gsub.replacement)) {
      if (length(gsub.pattern) != length(gsub.replacement))
        stop2("'gsub.replacement' and 'gsub.pattern' must have the same length")

      for (i in 1:length(gsub.pattern)) {
        posteriors_by_chain$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               posteriors_by_chain$variable)
        fit$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               fit$variable)
        coef_table2$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               coef_table2$variable)
        rownames(coef_table) <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               rownames(coef_table))
      }
    }

posteriors_by_chain$variable <-
      factor(posteriors_by_chain$variable,
             levels = sort(unique(posteriors_by_chain$variable), decreasing = TRUE))

    # define color for posterior distribution
    fill_values <-
      if (all(coef_table2$significance == "non-sig"))
        grDevices::adjustcolor(fill, alpha.f = 0.25) else
          if (all(coef_table2$significance == "sig"))
            grDevices::adjustcolor(fill, alpha.f = 0.5) else
              c(
                grDevices::adjustcolor(fill, alpha.f = 0.25),
                grDevices::adjustcolor(fill, alpha.f = 0.5)
              )

    if (!highlight)
      fill_values <- rep(grDevices::adjustcolor(fill, alpha.f = 0.5), length(fill_values))

    fit$significance <-
      sapply(fit$variable, function(x)
        coef_table2$significance[as.character(coef_table2$variable) == x])

    # choose effects to display
    if (!is.null(effects)){
      fit <- fit[grep(paste0("^", paste(effects, collapse = "$|^"), "$"), fit$variable), ]
      coef_table <- coef_table[grep(paste0("^", paste(effects, collapse = "$|^"), "$"), rownames(coef_table)), ]
      coef_table2 <- coef_table2[grep(paste0("^", paste(effects, collapse = "$|^"), "$"), coef_table2$variable), ]

      # order factors for plotting
      coef_table2 <- coef_table2[match(effects, coef_table2$variable), ]
      coef_table <- coef_table[match(effects, rownames(coef_table)), ]
      coef_table2$variable <- factor(coef_table2$variable, levels = effects)

      posteriors_by_chain <- posteriors_by_chain[grep(paste0("^", paste(effects, collapse = "$|^"), "$"), posteriors_by_chain$variable), ]
    }

    # define colors for point range in distribution plot
    if (highlight)
      col_pointrange <- ifelse(coef_table2$significance == "non-sig", "gray", "black") else col_pointrange <- rep("black", nrow(coef_table2))

    # order effects as in table
      fit$variable <- factor(fit$variable, levels = coef_table2$variable)

    # trick for getting own palette in ggplot2
      scale_color_discrete <- function(...) ggplot2::scale_color_manual(..., values= trace.palette(length(unique(posteriors_by_chain$chain))))

      on.exit(rm("scale_color_discrete"))

    # creat plots
    gg_distributions <-
      ggplot2::ggplot(data = fit, ggplot2::aes(y = variable, x = value, fill = significance)) +
      ggplot2::geom_vline(xintercept = 0,
                          col = "black",
                          lty = 2) +
      ggdist::stat_halfeye(
        ggplot2::aes(x = value),
        .width = c(.95),
        normalize = "panels",
        color = "transparent"
      ) +
      ggplot2::scale_fill_manual(values = fill_values, guide = 'none') +
      ggplot2::geom_point(data = coef_table2, col = col_pointrange) +
      ggplot2::geom_errorbar(data = coef_table2,
                             ggplot2::aes(xmin = `l-95% CI`, xmax = `u-95% CI`),
                             width = 0, col = col_pointrange) +
      ggplot2::facet_wrap(
        ~ variable,
        scales = "free_y",
        ncol = 1
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.ticks.length = ggplot2::unit(0, "pt"),
        plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme(
        panel.spacing.y = ggplot2::unit(0, "null")  # no vertical space between panels
      )


if (plot.area.prop != 1)
  gg_distributions <- gg_distributions + ggplot2::xlim(range(c(posteriors_by_chain$value, 0)) * plot.area.prop)

    if (trace){
    gg_traces <-
      ggplot2::ggplot(data = posteriors_by_chain, ggplot2::aes(x = iteration, y = value, color = chain)) +
      ggplot2::geom_line() +
     scale_color_discrete() +
      ggplot2::facet_wrap(
        ~ variable,
        scales = "free_y",
        ncol = 1,
        strip.position = "right"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      )

    gg <-
      cowplot::plot_grid(gg_distributions,
                         gg_traces,
                         ncol = 2,
                         rel_widths = c(1.8, 1))

    } else
      gg <- gg_distributions

    if (save){

      dir.create(file.path(dest.path, fit.name))


      if (height == "dynamic")
        height <- 3 + 0.4 * length(betas)
      if (height > 49) height <- 49

      cowplot::ggsave2(filename = file.path(dest.path, fit.name, "plot.jpeg"), plot = gg, width = width, height = height)

      saveRDS(object = list(fit_table = fit_table, coef_table = coef_table, graph = gg), file.path(dest.path, fit.name, "fit_table.RDS"))
      }
    if (return){
      output_list <- list(
        fit_table_html = html_format_fit_table(fit_table),
        coef_table = coef_table,
        coef_table_html = html_format_coef_table(coef_table, fill = fill,  highlight = highlight),
        plot = gg
        )
      return(output_list)
    } else {

      if (print.name)
        cat('\n\n## ', fit.name, '\n\n')

        # print fit summary table
        fit_table <- html_format_fit_table(fit_table)

        print(fit_table)

        # print estimates
        coef_table <- html_format_coef_table(coef_table, fill = fill,  highlight = highlight)

        # print fit result table
        print(coef_table)

         print(gg)
      }
  } else
    message("Folder already exists and overwrite = FALSE")
}
