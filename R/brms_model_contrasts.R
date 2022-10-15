brms_model_contrasts <- function(model, predictor, level.sep = " - ", xlab = "Effect size", n.subposts = 2000, gsub.pattern = NULL, gsub.replacement = NULL, fill = "#6DCD59FF"){

  # object for avoiding errors with ggplot functions when checking package
  significance <- unit <- value <- variable <- CI_high <- CI_low <- Hypothesis <- Parameter <- chain <- iteration <- margin <- NULL

  fit <- model$fit
  betas <- grep("^b_", names(fit@sim$samples[[1]]), value = TRUE)
  cats <- gsub("^b_", "", betas)
  contrsts <- paste(apply(utils::combn(cats, 2), 2, paste, collapse = " - "), "= 0")

  # fix  baseline level
  original_levels <- gsub(" ", "", paste0(predictor, unique(model$data[,predictor])))
  base_level <- setdiff(original_levels, cats)
  names(contrsts) <- gsub("Intercept", base_level, contrsts)
  contrsts <- gsub("- Intercept|Intercept - ", "", contrsts)

  names(contrsts) <- gsub(paste0(predictor,"|sensory_input| = 0"), "", names(contrsts))
  names(contrsts) <- gsub(" - ", level.sep, names(contrsts))

  if(!is.null(gsub.pattern) & !is.null(gsub.replacement)){

    if (length(gsub.pattern) != length(gsub.replacement))
      stop2("'gsub.replacement' and 'gsub.pattern' must have the same length")

    for(i in 1:length(gsub.pattern))
      names(contrsts) <- gsub(gsub.pattern[i], gsub.replacement[i], names(contrsts))
  }


  # evaluate hypothesis
  hyps <- brms::hypothesis(model, contrsts)

  hyp_table <- hyps$hypothesis[, c("Hypothesis", "Estimate", "Est.Error", "CI.Lower", "CI.Upper")]

  signif <- hyp_table[,"CI.Lower"] * hyp_table[,"CI.Upper"] > 0

  df1 <- kbl(hyp_table, row.names = TRUE, escape = FALSE, format = "html", digits = 3)

  df1 <- row_spec(kable_input = df1,row =  which(hyp_table$CI.Lower * hyp_table$CI.Upper > 0), background = grDevices::adjustcolor(fill, alpha.f = 0.3))

  kbl1 <- kable_styling(df1, bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE, font_size = 12)

  # subsample posteriors
  xdrws <- brms::as_draws(hyps$samples)

  names(xdrws)[1:length(contrsts)] <- names(contrsts)

  # only apply thinning if length of posterior < n.subposts
  xdrws <- posterior::subset_draws(x = xdrws, variable = names(contrsts))

  if (round(length(xdrws[[1]][[1]]) / n.subposts, 0) >= 2)
    xdrws <- posterior::thin_draws(xdrws, round(length(xdrws[[1]][[1]]) / n.subposts, 0))

  merged_xdrws <- posterior::merge_chains(xdrws)
  sub_posts <- as.data.frame(merged_xdrws)[, names(contrsts)]
  names(sub_posts) <- names(contrsts)

  hdis <- t(sapply(names(contrsts), function(y)   HDInterval::hdi(sub_posts[, colnames(sub_posts) == y]))
  )

  hyp_table$CI_low <- round(hyp_table$CI.Lower, digits = 3)
  hyp_table$CI_high <- round(hyp_table$CI.Upper, digits = 3)
  hyp_table$CI.Lower <- hyp_table$CI.Upper <- NULL

  out <- lapply(names(contrsts), function(y)  data.frame(variable = y, value = sort(sub_posts[, colnames(sub_posts) == y], decreasing = FALSE)))

  posteriors <- do.call(rbind, out)
  posteriors$Hypothesis <- factor(posteriors$variable, levels = sort(unique(posteriors$variable)))

  coef_table2 <- hyp_table
  coef_table2$value <- coef_table2$Estimate
  coef_table2$significance <- ifelse(coef_table2$CI_low * coef_table2$CI_high > 0, "sig", "non-sig")
  coef_table2$significance <- factor(coef_table2$significance, levels = c("non-sig", "sig"))

  col_pointrange <- if(all(coef_table2$significance == "non-sig")) "gray" else
    if(all(coef_table2$significance == "sig")) "black" else
      c("gray", "black")

  posteriors$significance <- sapply(posteriors$Hypothesis, function(x) coef_table2$significance[coef_table2$Hypothesis == x])

  gg_distribution <- ggplot2::ggplot(data = posteriors, ggplot2::aes(y = Hypothesis, x = value, fill = significance)) +
    ggplot2::geom_vline(xintercept = 0, col = "black", lty = 2) +
    ggdist::stat_halfeye(ggplot2::aes(x = value), .width = c(.95),  normalize = "panels", color = "transparent") +
    ggplot2::scale_fill_manual(values = c(grDevices::adjustcolor(fill, alpha.f = 0.25), grDevices::adjustcolor(fill, alpha.f = 0.5))) +
    ggplot2::geom_point(data = coef_table2) +
    ggplot2::geom_errorbar(data = coef_table2, ggplot2::aes(xmin = CI_low, xmax = CI_high), width = 0) +
    ggplot2::scale_color_manual(values = col_pointrange) +
    ggplot2::theme(axis.ticks.length = unit(0, "pt"), plot.margin = margin(0,0,0,0,"pt"), legend.position="none",
                   strip.background = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank()) +
    ggplot2::labs(x = xlab, y = "Contrasts") +
    ggplot2::theme_classic()

  print(kbl1)
  print(gg_distribution)

}
