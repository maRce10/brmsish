
# reads and prints summaries saved by html_summary()
read_summary <- function(path = ".", fill = "#6DCD59FF"){

    # Read model output
    mod <- readRDS(file.path(i, "model_table.RDS"))

    # print model formula
    # cat('\n\n')
    cat('\n\n## ', basename(path), '\n\n')

    # print model summary table
    model_table <- html_format_model_table(mod$model_table)

    print(model_table)

    # print estimates
    coef_table <- html_format_coef_table(mod$coef_table, fill = fill)

    # print model result table
    print(coef_table)

    # plot
    path <- normalizePath(path)
    cat("![](", file.path(path, "plot.jpeg"), ")")
}


html_format_model_table <- function(x){

  x$diverg_transitions <-
    ifelse(
      x$diverg_transitions > 0,
      kableExtra::cell_spec(
        x$diverg_transitions,
        "html",
        color = "white",
        background = "red",
        bold = TRUE,
        font_size = 12
      ),
      kableExtra::cell_spec(x$diverg_transitions, "html")
    )

  x$`rhats > 1.05` <-
    ifelse(
      x$rhats...1.05 > 0,
      kableExtra::cell_spec(
        x$rhats...1.05,
        "html",
        color = "white",
        background = "red",
        bold = TRUE,
        font_size = 12
      ),
      kableExtra::cell_spec(x$rhats...1.05, "html")
    )

  # remove all rhats column
  x$rhats...1.05 <- NULL

  # sort columns
  x <- x[, c("priors",	"formula",	"iterations", "chains", "thinning", "warmup",	"diverg_transitions", "rhats > 1.05",	"min_bulk_ESS", "min_tail_ESS", "seed")]


  x_kbl <-
    kableExtra::kbl(
      x,
      row.names = TRUE,
      escape = FALSE,
      format = "html",
      digits = 3
    )

  x_kbl <-
    kableExtra::kable_styling(
      x_kbl,
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      font_size = 12
    )

}


html_format_coef_table <- function(x, fill){


  x$Rhat <- round(x$Rhat, 3)

  x$Rhat <-
    ifelse(
      x$Rhat > 1.05,
      kableExtra::cell_spec(
        x$Rhat,
        "html",
        color = "white",
        background = "red",
        bold = TRUE,
        font_size = 12
      ),
      kableExtra::cell_spec(x$Rhat, "html")
    )

  x_kbl <-
    kableExtra::kbl(
      x,
      row.names = TRUE,
      escape = FALSE,
      format = "html",
      digits = 3
    )

  x_kbl <-
    kableExtra::row_spec(
      kable_input = x_kbl,
      row =  which(x$`l-95% CI` * x$`u-95% CI` > 0),
      background = grDevices::adjustcolor(fill, alpha.f = 0.3)
    )

  x_kbl <-
    kableExtra::kable_styling(
      x_kbl,
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      font_size = 12
    )
}
