# prints a summary of brm model results. It includes a model table, a coefficient table and a posterior distribution halfeye graph next to a chain trace plot. Tables are produce in html format so the output is nicely printed when knit in Rmarkdown files.

# remember to add 'results = 'as.is' and warning = FALSE to chunk options and adjust fig.width
html_summary <-
  function(model = NULL,
           gsub.pattern = NULL, # a vector with character strings to be replaced
           gsub.replacement = NULL,  # a vector with character strings to use for replacement
           xlab = "Effect size",
           n.posterior = 2000, # number of posterior samples to use for plotting
           model.name = NULL, # for adding it as a title
           read.file = NULL, # name of file to be read, if supplied 'model' is ignored
           plot.area.prop = 1, # proportion of the total effect size range (given by the minimum and maximum of all posterior samples) to plot (sort of xlim)
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
           print.name = TRUE
           ) {

     # object for avoiding errors with ggplot functions when checking package
    significance <-
      value <-
      variable <-
      CI_high <-
      CI_low <-
      Hypothesis <- Parameter <- chain <- iteration <- NULL

    if (is.null(model) & is.null(read.file))
      stop("either 'model' or 'read.file' must be supplied")

    if (!is.null(model) & is.null(model.name))
      model.name <- deparse(substitute(model)) else
        if (!is.null(read.file) & is.null(model.name))
          model.name <- gsub("\\.rds$", "", basename(read.file), ignore.case = TRUE)

    # skip everything if save TRUE and folder exists
    if (!dir.exists(file.path(dest.path, model.name)) & save | !save | save & overwrite){

    if (is.null(model) & !is.null(read.file))
      model <- readRDS(read.file)
    variables <- posterior::variables(model)
    incl_classes <- c(
      "b", "bs", "bcs", "bsp", "bmo", "bme", "bmi", "bm",
      brms:::valid_dpars(model), "delta", "lncor", "rescor", "ar", "ma", "sderr",
      "cosy", "cortime", "lagsar", "errorsar", "car", "sdcar", "rhocar",
      "sd", "cor", "df", "sds", "sdgp", "lscale", "simo"
    )
    incl_regex <- paste0("^", brms:::regex_or(incl_classes), "(_|$|\\[)")
    variables <- variables[grepl(incl_regex, variables)]
    betas <- grep("^b_", variables, value = TRUE)

    # remove intercept betas
    if (remove.intercepts)
      betas <- grep("b_Intercept", betas, value = TRUE, invert = TRUE)

    iterations <- model$fit@sim$iter
    chains <- posterior::nchains(model)
    warmup <- model$fit@sim$warmup
    mod_formula <- as.character(model$formula[1])
    diverg_transitions <- sum(brms::nuts_params(model, pars = "divergent__")$Value)
    priors <- paste(apply(brms::prior_summary(model, all = FALSE)[, 2:1], 1, paste, collapse = "-"), collapse = "\n")
    seed <- model$fit@stan_args[[1]]$seed
    thinning <- model$fit@stan_args[[1]]$thin

    # replace model with draws (to avoid having several huge objects)
    model <- posterior::as_draws_array(model, variable = betas)

    coef_table <- draw_summary(model, variables = betas, probs = c(0.025, 0.975), robust = robust)

    model_table <-
      data.frame(
        priors = priors,
        formula = mod_formula,
        iterations = iterations,
        chains = chains,
        thinning = thinning,
        warmup = warmup,
        diverg_transitions = diverg_transitions,
        `rhats > 1.05` = sum(coef_table$Rhat > 1.05),
        min_bulk_ESS = min(coef_table$Bulk_ESS),
        min_tail_ESS = min(coef_table$Tail_ESS),
        seed = seed
      )

    # thin before getting chain data for trace plot
    if (round(model_table$iterations / n.posterior, 0) >= 2)
      model <-
      posterior::thin_draws(model, model_table$iterations / n.posterior, 0)

    # data by chain for trace plot
    posteriors_by_chain <- do.call(rbind, lapply(betas, function(x) {

      X <- as.data.frame(model[,, x])
      names(X) <- paste("chain", 1:ncol(X))
      X <- stack(X)
      X$variable <- x
      X$iteration <-
        round(seq(1, model_table$iterations, length.out = nrow(X) / model_table$chains))
      names(X) <- c("value", "chain", "variable", "iteration")
      return(X)
    }))

    # merge chains
    model <- as.data.frame(posterior::merge_chains(model))
    names(model) <- betas

    model <- do.call(rbind, lapply(betas, function(y)
      data.frame(
        variable = y,
        value = sort(model[, colnames(model) == y], decreasing = FALSE)
      )))

    model$variable <-
      factor(model$variable, levels = sort(unique(model$variable)))

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
        model$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               model$variable)
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

# define colors for point range in distribution plot
if (highlight)
col_pointrange <- ifelse(coef_table2$significance == "non-sig", "gray", "black")
else col_pointrange <- rep("black", nrow(coef_table2))

    # define color for posterior distribution
    # fill_values <- rep(grDevices::adjustcolor(fill, alpha.f = 0.5), nrow(coef_table2))
    #
    # if (highlight)
    # fill_values <- ifelse(coef_table2$significance == "no-sig", grDevices::adjustcolor(fill, alpha.f = 0.25), fill_values)

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

    model$significance <-
      sapply(model$variable, function(x)
        coef_table2$significance[as.character(coef_table2$variable) == x])

    # choose effects to display
    if (!is.null(effects)){
      model <- model[grep(paste(effects, collapse = "|"), model$variable), ]
      coef_table <- coef_table[grep(paste(effects, collapse = "|"), rownames(coef_table)), ]
      coef_table2 <- coef_table2[grep(paste(effects, collapse = "|"), coef_table2$variable), ]
      posteriors_by_chain <- posteriors_by_chain[grep(paste(effects, collapse = "|"), posteriors_by_chain$variable), ]
    }

    # order effects as in table
      model$variable <- factor(model$variable, levels = coef_table2$variable)

    # trick for getting own palette in ggplot2
      scale_color_discrete <- function(...) ggplot2::scale_color_manual(..., values= trace.palette(length(unique(posteriors_by_chain$chain))))

      on.exit(rm("scale_color_discrete"))

    # creat plots
    gg_distributions <-
      ggplot2::ggplot(data = model, ggplot2::aes(y = variable, x = value, fill = significance)) +
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
      ggplot2::labs(x = "Effect size", y = "Parameter")


if (plot.area.prop != 1)
  gg_distributions <- gg_distributions + ggplot2::xlim(range(c(posteriors_by_chain$value, 0)) * plot.area.prop)

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

    if (save){

      dir.create(file.path(dest.path, model.name))


      if (height == "dynamic")
        height <- 3 + 0.4 * length(betas)
      if (height > 49) height <- 49

      cowplot::ggsave2(filename = file.path(dest.path, model.name, "plot.jpeg"), plot = gg, width = width, height = height)
    }

    # save output
    if (save)
      saveRDS(object = list(model_table = model_table, coef_table = coef_table, graph = gg), file.path(dest.path, model.name, "model_table.RDS")) else {

      if (print.name)
        cat('\n\n## ', model.name, '\n\n')

        # print model summary table
        model_table <- html_format_model_table(model_table)

        print(model_table)

        # print estimates
        coef_table <- html_format_coef_table(coef_table, fill = fill,  highlight = highlight)

        # print model result table
        print(coef_table)

         print(gg)
      }
  } else
    message("Folder already exists and overwrite = FALSE")
}
