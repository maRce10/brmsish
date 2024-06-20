
draw_extended_summary <- function(draws,
                                  name = NULL,
                                  highlight = TRUE,
                                  fill = "#6DCD59FF",
                                  remove.intercepts = FALSE,
                                  by = NULL,
                                  gsub.pattern = NULL,
                                  gsub.replacement = NULL,
                                  xlab = "Effect size",
                                  ylab = "Parameter") {

  # create objects just to avoid errors with ggplot functions when checking package
  position_dodge <- org.variable <- value <- significance <- `l-95% CI` <- `u-95% CI` <- theme <- unit <- NULL

  if (!is.null(by)) {
    levels <- draws[, by]

    unique_levels <- unique(levels)
  } else
    unique_levels <- ".A"

  # run loop over each level (or just the single "level")
  results_list <- lapply(unique_levels, function(x) {
    if (identical(unique_levels, ".A")) {
      # keep only betas
      draws <- draws[, grep("^b_", names(draws), value = TRUE)]
    } else {
      draws <-
        draws[draws[, by, drop = TRUE] == x, grep("^b_", names(draws), value = TRUE)]
    }

    # remove intercept betas
    if (remove.intercepts)
      draws <-
        draws[, grep("^b_Intercept",
                     names(draws),
                     value = TRUE,
                     invert = TRUE), drop = FALSE]

    # compute model-averaged posteriors of overlapping parameters
    coef_table <-
      posterior::summarise_draws(draws,
                                 median,
                                 ~ quantile(.x, probs = c(0.025, 0.975)),
                                 posterior::default_convergence_measures())

    names(coef_table)[3:4] <- c("l-95% CI", "u-95% CI")

    coef_table$value <- coef_table$median
    coef_table$significance <-
      ifelse(coef_table$`l-95% CI` * coef_table$`u-95% CI` > 0,
             "sig",
             "non-sig")
    coef_table$significance <-
      factor(coef_table$significance, levels = c("non-sig", "sig"))


    if (ncol(draws) > 1)
      sdraws <-
      stack(draws[, grep("^b_", names(draws), value = TRUE)], ind = "levels")
    else
      sdraws <-
      data.frame(value = draws[, 1], variable = names(draws))

    names(sdraws) <- c("value", "variable")

    sdraws$significance <-
      sapply(sdraws$variable, function(x)
        coef_table$significance[as.character(coef_table$variable) == x][1])

    # add level
    # if (!identical(unique_levels, ".A")) {
    # to coef_table
    coef_table$.new.column <- x
    names(coef_table)[ncol(coef_table)] <- "levels"

    #to sdraws
    sdraws$.new.column <- x
    names(sdraws)[ncol(sdraws)] <- "levels"
    # }

    output <- list(coef_table = coef_table, sdraws = sdraws)

    return(output)
  })

  # put both results into a single data frame
  coef_table <- do.call(rbind, lapply(results_list, "[[", 1))
  sdraws <- do.call(rbind, lapply(results_list, "[[", 2))

  if (!is.null(gsub.pattern) & !is.null(gsub.replacement)) {
    if (length(gsub.pattern) != length(gsub.replacement))
      stop2("'gsub.replacement' and 'gsub.pattern' must have the same length")

    for (i in 1:length(gsub.pattern)) {
      sdraws$variable <-
        gsub(pattern = gsub.pattern[i],
             replacement = gsub.replacement[i],
             sdraws$variable)
      coef_table$variable <-
        gsub(pattern = gsub.pattern[i],
             replacement = gsub.replacement[i],
             coef_table$variable)
    }
  }


  # duplicate variable column
  coef_table <- as.data.frame(coef_table)
  coef_table$org.variable <- coef_table$variable
  sdraws$org.variable <- sdraws$variable

  if (!identical(unique_levels, ".A")){
    coef_table$variable <- paste(coef_table$levels, coef_table$variable, sep = "-")
    coef_table <- coef_table[order(coef_table$org.variable), ]

    sdraws$variable <- paste(sdraws$levels, sdraws$variable, sep = "-")
  }

  rownames(coef_table) <- coef_table$variable

  fill_df <- data.frame(level = unique_levels, fill = fill)

  if (!identical(unique_levels, ".A"))
    coef_table$fill_values <- sapply(coef_table$levels, function(x) fill_df$fill[fill_df$level == x]) else
      coef_table$fill_values <- fill[1]

  if (highlight)
  {
    coef_table$col_pointrange <-
      ifelse(coef_table$significance == "non-sig", "gray", "black")

    # coef_table$fill_values <- ifelse(coef_table$significance == "no-sig", grDevices::adjustcolor(col = coef_table$fill_values, alpha.f = 0.5), coef_table$fill_values)
    sdraws$significance <- sapply(sdraws$variable, function(x) coef_table$significance[coef_table$variable == x])
  }  else  {
    coef_table$col_pointrange <- rep("black", nrow(coef_table))
    sdraws$significance <- 1
  }

  pd <- position_dodge(width = 0.05)

  # creat plots
  gg_distributions <-
    ggplot2::ggplot(data = sdraws, ggplot2::aes(y = org.variable, x = value, fill = levels, alpha = if(highlight) significance else NULL)) +
    ggplot2::geom_vline(xintercept = 0,
                        col = "black",
                        lty = 2) +
    ggdist::stat_halfeye(
      ggplot2::aes(x = value),
      .width = c(.95),
      normalize = "panels",
      color = "transparent", position = pd) +
    ggplot2::scale_alpha_manual(values = c(0.4, 0.8), guide = 'none') +
    ggplot2::scale_fill_manual(values = if (!identical(unique_levels, ".A")) as.vector(coef_table$fill_values) else fill[1]) +
    ggplot2::geom_point(data = coef_table, position = pd) +
    ggplot2::geom_errorbar(
      data = coef_table,
      ggplot2::aes(xmin = `l-95% CI`, xmax = `u-95% CI`),
      width = 0,
      position = pd
    ) +
    ggplot2::facet_wrap(~ org.variable,
                        scales = "free_y",
                        ncol = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.length = ggplot2::unit(0, "pt"),
      plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, fill = "Effect") +
    theme(
      panel.spacing.y = unit(0, "null")  # no vertical space between panels
    )

    coef_table$variable <-
      coef_table$significance <-  coef_table$value <- NULL
    names(coef_table) <-
      c("Estimate",
        "l-95% CI",
        "u-95% CI",
        "Rhat",
        "Bulk_ESS",
        "Tail_ESS")

    coef_table <- coef_table[, c("Estimate",
                                 "l-95% CI",
                                 "u-95% CI",
                                 "Rhat",
                                 "Bulk_ESS",
                                 "Tail_ESS")]

  if (!is.null(name))
    cat('\n\n## ', name, '\n\n')

  html_coef_table <-
    html_format_coef_table(coef_table, fill = fill,  highlight = highlight)

  print(html_coef_table)

  print(gg_distributions)
}
