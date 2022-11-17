# prints a summary of brm model results. It includes a model table, a coefficient table and a posterior distribution halfeye graph next to a chain trace plot. Tables are produce in html format so the output is nicely printed when knit in Rmarkdown files.

# remember to add 'results = 'as.is' and warning = FALSE to chunk options and adjust fig.width
html_summary <-
  function(model = NULL,
           gsub.pattern = NULL, # a vector with character strings to be replaced
           gsub.replacement = NULL,  # a vector with character strings to use for replacement
           xlab = "Effect size",
           n.posterior = 2000, # number of posterior samples to use for plotting
           model_name = NULL, # for adding it as a title
           read.file = NULL, # name of file to be read, if supplied 'model' is ignored
           plot.area.prop = 1, # proportion of the total effect size range (given by the minimum and maximum of all posterior samples) to plot (sort of xlim)
           remove.intercepts = FALSE,
           fill = "#6DCD59FF",
           trace.palette = viridis,
           effects = NULL) {

     # object for avoiding errors with ggplot functions when checking package
    significance <-
      value <-
      variable <-
      CI_high <-
      CI_low <-
      Hypothesis <- Parameter <- chain <- iteration <- NULL

    if (is.null(model) & !is.null(read.file))
      model <- readRDS(read.file)

    # extract info from model
    summ <- summary(model)$fixed

    if (remove.intercepts)
      summ <-
      summ[grep("^Intercept", rownames(summ), invert = TRUE),]


    fit <- model$fit
    betas <- grep("^b_", names(fit@sim$samples[[1]]), value = TRUE)

    # remove intercept betas
    if (remove.intercepts)
      betas <- grep("b_Intercept", betas, value = TRUE, invert = TRUE)

    # subsample posteriors
    xdrws <- brms::as_draws(model)

    # only apply thinning if length of posterior < n.posterior
    if (round(length(xdrws[[1]][[1]]) / n.posterior, 0) >= 2)
      xdrws <-
      posterior::thin_draws(xdrws, round(length(xdrws[[1]][[1]]) / n.posterior, 0))

    xdrws <- posterior::subset_draws(x = xdrws, variable = betas)
    sub_posts_by_chain_list <- lapply(1:length(xdrws), function(x) {
      X <- as.data.frame(xdrws[[x]])
      X$chain <- paste("chain", x)
      return(X)
    })

    sub_posts_by_chain <- do.call(rbind, sub_posts_by_chain_list)

    merged_xdrws <- posterior::merge_chains(xdrws)
    sub_posts <- as.data.frame(merged_xdrws)
    names(sub_posts) <- betas

    coef_table <- data.frame(summ)
    coef_table <-
      coef_table[, c("Estimate", "Rhat", "Bulk_ESS", "Tail_ESS", "l.95..CI", "u.95..CI")]

    # add priors to model table
    pt <- prior_summary(model)
    b_prior <- pt$prior[pt$class == "b" & pt$coef == ""]
    b_prior <- if (b_prior == "")
      "flat" else
      b_prior

    sd_prior <- unique(pt$prior[pt$class == "sd" & pt$coef == ""])
    sd_prior <- if (length(sd_prior) > 1)
      sd_prior[sd_prior != ""] else
      "flat"
    sd_prior <- if (sd_prior == "")
      "flat" else
      sd_prior

    model_table <-
      data.frame(
        b_prior,
        sd_prior,
        iterations = fit@stan_args[[1]]$iter,
        chains = length(attr(fit, "stan_args")),
        thinning = fit@stan_args[[1]]$thin,
        warmup = fit@stan_args[[1]]$warmup
      )

    np <- brms::nuts_params(model)
    model_table$diverg_transitions <-
      sum(subset(np, Parameter == "divergent__")$Value)
    model_table$`rhats > 1.05` <-
      sum(stats::na.omit(brms::rhat(model)) > 1.05)

    model_table$min_bulk_ESS <- min(coef_table$Bulk_ESS)
    model_table$min_tail_ESS <- min(coef_table$Tail_ESS)

    model_table$seed <- fit@stan_args[[1]]$seed

    coef_table <- as.data.frame(coef_table)
    coef_table$Rhat <- round(coef_table$Rhat, digits = 3)
    coef_table$CI_low <-
      round(unlist(coef_table$l.95..CI), digits = 3)
    coef_table$CI_high <-
      round(unlist(coef_table$u.95..CI), digits = 3)
    coef_table$l.95..CI <- coef_table$u.95..CI <- NULL

    out <-
      lapply(betas, function(y)
        data.frame(
          variable = y,
          value = sort(sub_posts[, colnames(sub_posts) == y], decreasing = FALSE)
        ))

    posteriors <- do.call(rbind, out)
    posteriors$variable <-
      factor(posteriors$variable, levels = sort(unique(posteriors$variable)))

    names(sub_posts_by_chain)[1:(ncol(sub_posts_by_chain) - 1)] <-
      betas

    out2 <- lapply(betas, function(y)  {
      X <-
        data.frame(variable = y,
                   chain = sub_posts_by_chain[, "chain"],
                   value = sub_posts_by_chain[, y])
      X$iteration <-
        round(seq(1, fit@stan_args[[1]]$iter, length.out = nrow(X) / length(attr(
          fit, "stan_args"
        ))))
      return(X)
    })

    posteriors_by_chain <- do.call(rbind, out2)

    coef_table2 <- coef_table
    coef_table2$variable <-
      factor(paste0("b_", rownames(coef_table2)))
    coef_table2$value <- coef_table2$Estimate
    coef_table2$significance <-
      ifelse(coef_table2$CI_low * coef_table2$CI_high > 0, "sig", "non-sig")
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
        posteriors$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               posteriors$variable)
        coef_table2$variable <-
          gsub(pattern = gsub.pattern[i],
               replacement = gsub.replacement[i],
               coef_table2$variable)

      }
    }


    if (!is.null(gsub.pattern) & !is.null(gsub.replacement)) {

    }

    posteriors_by_chain$variable <-
      factor(posteriors_by_chain$variable,
             levels = sort(unique(posteriors_by_chain$variable), decreasing = TRUE))

    col_pointrange <-
      if (all(coef_table2$significance == "non-sig"))
        "gray" else
      if (all(coef_table2$significance == "sig"))
        "black" else
      c("gray", "black")

    fill_values <-
      if (all(coef_table2$significance == "non-sig"))
        grDevices::adjustcolor(fill, alpha.f = 0.25) else
          if (all(coef_table2$significance == "sig"))
            grDevices::adjustcolor(fill, alpha.f = 0.5) else
              c(
                grDevices::adjustcolor(fill, alpha.f = 0.25),
                grDevices::adjustcolor(fill, alpha.f = 0.5)
              )


    posteriors$significance <-
      sapply(posteriors$variable, function(x)
        coef_table2$significance[as.character(coef_table2$variable) == x])

    # choose effects to display
    if (!is.null(effects)){
      posteriors <- posteriors[grep(paste(effects, collapse = "|"), posteriors$variable), ]
      coef_table <- coef_table[grep(paste(effects, collapse = "|"), rownames(coef_table)), ]
      coef_table2 <- coef_table2[grep(paste(effects, collapse = "|"), coef_table2$variable), ]
      posteriors_by_chain <- posteriors_by_chain[grep(paste(effects, collapse = "|"), posteriors_by_chain$variable), ]
    }

    # order effects as in table
      posteriors$variable <- factor(posteriors$variable, levels = coef_table2$variable)

    # trick for getting own palette in ggplot2
      scale_color_discrete <- function(...) scale_color_manual(..., values= trace.palette(length(unique(posteriors_by_chain$chain))))

      on.exit(rm("scale_color_discrete"))


    # creat plots
    gg_dists <-
      ggplot2::ggplot(data = posteriors, ggplot2::aes(y = variable, x = value, fill = significance)) +
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
      ggplot2::geom_point(data = coef_table2) +
      ggplot2::geom_errorbar(data = coef_table2,
                             ggplot2::aes(xmin = CI_low, xmax = CI_high),
                             width = 0) +
      ggplot2::scale_color_manual(values = col_pointrange) +
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
      ggplot2::labs(x = "Effect size", y = "Parameter") +

      ggplot2::xlim(range(c(posteriors_by_chain$value, 0)) * plot.area.prop)

    gg_traces <-
      ggplot2::ggplot(data = posteriors_by_chain, ggplot2::aes(x = iteration, y = value, color = chain)) +
      ggplot2::geom_line() +
      # ggplot2::scale_color_viridis_d(alpha = 0.7,
      #                                begin = 0.2,
      #                                end = 0.9) +
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
      cowplot::plot_grid(gg_dists,
                         gg_traces,
                         ncol = 2,
                         rel_widths = c(2, 1))

    if (!is.null(gsub.pattern) & !is.null(gsub.replacement))
      rownames(coef_table) <-
      gsub(pattern = gsub.pattern,
           replacement = gsub.replacement,
           rownames(coef_table))

    coef_table$Rhat <-
      ifelse(
        coef_table$Rhat > 1.05,
        cell_spec(
          coef_table$Rhat,
          "html",
          color = "white",
          background = "red",
          bold = TRUE,
          font_size = 12
        ),
        cell_spec(coef_table$Rhat, "html")
      )

    signif <- coef_table[, "CI_low"] * coef_table[, "CI_high"] > 0

    model_table$diverg_transitions <-
      ifelse(
        model_table$diverg_transitions > 0,
        cell_spec(
          model_table$diverg_transitions,
          "html",
          color = "white",
          background = "red",
          bold = TRUE,
          font_size = 12
        ),
        cell_spec(model_table$diverg_transitions, "html")
      )

    model_table$`rhats > 1.05` <-
      ifelse(
        model_table$`rhats > 1.05` > 0,
        cell_spec(
          model_table$`rhats > 1.05`,
          "html",
          color = "white",
          background = "red",
          bold = TRUE,
          font_size = 12
        ),
        cell_spec(model_table$`rhats > 1.05`, "html")
      )


    df1 <-
      kbl(
        model_table,
        row.names = TRUE,
        escape = FALSE,
        format = "html",
        digits = 3
      )

    df1 <-
      kable_styling(
        df1,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE,
        font_size = 12
      )


    df2 <-
      kbl(
        coef_table,
        row.names = TRUE,
        escape = FALSE,
        format = "html",
        digits = 3
      )

    df2 <-
      row_spec(
        kable_input = df2,
        row =  which(coef_table$CI_low * coef_table$CI_high > 0),
        background = grDevices::adjustcolor(fill, alpha.f = 0.3)
      )

    df2 <-
      kable_styling(
        df2,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE,
        font_size = 12
      )

    if (!is.null(model_name))
      cat(paste('<font size="4"><b>', model_name, '</b></font><br>'))

    cat(paste('<font size="3"><b>', model$formula[1], '</b></font>'))

    print(df1)
    print(df2)

    print(gg)
  }
