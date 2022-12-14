# Estimates hypothesis testing for all pairwise comparisons of levels from a categorical predictor. The function brms::hypothesis is used internally.

contrasts <-
  function(model, # a brms model
           predictor, # name of categorical predictor
           level.sep = " - ", # string to use for separating level names
           xlab = "Effect size",
           gsub.pattern = NULL, # a vector with character strings to be replaced
           gsub.replacement = NULL,  # a vector with character strings to use for replacement
           n.posterior = 2000, # number of posterior samples to use for plotting
           fill = "#6DCD59FF",
           sort.levels = NULL,
           html.table = FALSE,
           plot = FALSE
           ) {

    # object for avoiding errors with ggplot functions when checking package
    significance <-
      value <-
      variable <-
      CI_high <-
      CI_low <-
      Hypothesis <- Parameter <- chain <- iteration <- NULL

    fit <- model$fit
    betas <- grep("^b_", names(fit@sim$samples[[1]]), value = TRUE)
    model_levels <- gsub("^b_", "", betas)
    model_levels <- grep(paste0("^Intercept$|", predictor), model_levels, value = TRUE)

    # fix  baseline level
    original_levels <-
      gsub(" |&", "", paste0(predictor, unique(model$data[, predictor])))
    base_level <- setdiff(original_levels, model_levels)

    # get levels
    pred_levels <- as.character(unique(model$data[, predictor]))

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
    hyps <- brms::hypothesis(model, contrsts)

    hyp_table <-
      hyps$hypothesis[, c("Hypothesis",
                          "Estimate",
                          "Est.Error",
                          "CI.Lower",
                          "CI.Upper")]

    hyp_table$Estimate <- hyp_table$Estimate * levels_df$sign
    hyp_table$CI.Lower <- hyp_table$CI.Lower * levels_df$sign
    hyp_table$CI.Upper <- hyp_table$CI.Upper * levels_df$sign

    if (html.table){
      signif <- hyp_table[, "CI.Lower"] * hyp_table[, "CI.Upper"] > 0

      html_table <-
      kbl(
        hyp_table,
        row.names = TRUE,
        escape = FALSE,
        format = "html",
        digits = 3
      )

    html_table <-
      row_spec(
        kable_input = html_table,
        row =  which(hyp_table$CI.Lower * hyp_table$CI.Upper > 0),
        background = grDevices::adjustcolor(fill, alpha.f = 0.3)
      )

    html_table <-
      kable_styling(
        html_table,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width = FALSE,
        font_size = 12
      )

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

    hyp_table$CI_low <- round(hyp_table$CI.Lower, digits = 3)
    hyp_table$CI_high <- round(hyp_table$CI.Upper, digits = 3)
    hyp_table$CI.Lower <- hyp_table$CI.Upper <- NULL

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
      ggplot2::scale_color_manual(values = col_pointrange) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.ticks.length = ggplot2::unit(0, "pt"),
        plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
        legend.position = "none",
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = xlab, y = "Contrasts")


    print(gg_distribution)
}
    if (!html.table)
    return(hyp_table)

  }
