
# reads and prints summaries saved by html_summary()
read_summary <- function(path = ".", fill = "#6DCD59FF", highlight = FALSE){

    # Read model output
    mod <- readRDS(file.path(i, "model_table.RDS"))

    # print model formula
    # cat('\n\n')
    cat('\n\n## ', basename(path), '\n\n')

    # print model summary table
    model_table <- html_format_model_table(mod$model_table)

    print(model_table)

    # print estimates
    coef_table <- html_format_coef_table(mod$coef_table, fill = fill, highlight = highlight)

    # print model result table
    print(coef_table)

    # plot
    path <- normalizePath(path)
    cat("![](", file.path(path, "plot.jpeg"), ")")
}
