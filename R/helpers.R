# internal functions not to be called by users
# stop function that doesn't print call
stop2 <- function (...)
{
  stop(..., call. = FALSE)
}

# internal function, not to be called by users. It is a modified version of pbapply::pblapply
# that allows to define internally if progress bar would be used (pbapply::pblapply uses pboptions to do this)
#last modification on oct-31-2022 (MAS)
#'
pblapply_brmsish_int <- function(X, FUN, cl = 1, pbar = TRUE, ...) {

  # conver parallel 1 to null
  if (!inherits(cl, "cluster"))
    if (cl == 1) cl <- NULL

  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X))
    X <- as.list(X)
  if (!length(X))
    return(lapply(X, FUN, ...))
  if (!is.null(cl)) {
    if (.Platform$OS.type == "windows") {
      if (!inherits(cl, "cluster"))
        cl <- NULL
    } else {
      if (inherits(cl, "cluster")) {
        if (length(cl) < 2L)
          cl <- NULL
      } else {
        if (cl < 2)
          cl <- NULL
      }
    }
  }

  if (is.null(cl)) {
    if (!pbar)
      return(lapply(X, FUN, ...))
    Split <- pbapply::splitpb(length(X), 1L, nout = 100)
    B <- length(Split)
    pb <- pbapply::startpb(0, B)
    on.exit(pbapply::closepb(pb), add = TRUE)
    rval <- vector("list", B)
    for (i in seq_len(B)) {
      rval[i] <- list(lapply(X[Split[[i]]], FUN, ...))
      pbapply::setpb(pb, i)
    }
  } else {
    if (inherits(cl, "cluster")) {
      PAR_FUN <- parallel::parLapply
      if (pbar)
        return(PAR_FUN(cl, X, FUN, ...))
      Split <- pbapply::splitpb(length(X), length(cl), nout = 100)
      B <- length(Split)
      pb <- pbapply::startpb(0, B)
      on.exit(pbapply::closepb(pb), add = TRUE)
      rval <- vector("list", B)
      for (i in seq_len(B)) {
        rval[i] <- list(PAR_FUN(cl, X[Split[[i]]], FUN,
                                ...))
        pbapply::setpb(pb, i)
      }
    } else {
      if (!pbar)
        return(parallel::mclapply(X, FUN, ..., mc.cores = as.integer(cl)))
      Split <- pbapply::splitpb(length(X), as.integer(cl), nout = 100)
      B <- length(Split)
      pb <- pbapply::startpb(0, B)
      on.exit(pbapply::closepb(pb), add = TRUE)
      rval <- vector("list", B)
      for (i in seq_len(B)) {
        rval[i] <- list(parallel::mclapply(X[Split[[i]]],
                                           FUN, ..., mc.cores = as.integer(cl)))
        pbapply::setpb(pb, i)
      }
    }
  }
  rval <- do.call(c, rval, quote = TRUE)
  names(rval) <- names(X)
  rval
}

# prints tables with fit diagnostics in html format
html_format_fit_table <- function(x){

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

  # relabel rhats column
  if (!is.null(x$rhats...1.05))
    x$`rhats > 1.05` <- x$rhats...1.05

  # remove rhats column
  x$rhats...1.05 <- NULL

  x$`rhats > 1.05` <-
    ifelse(
      x$`rhats > 1.05` > 0,
      kableExtra::cell_spec(
        x$`rhats > 1.05`,
        "html",
        color = "white",
        background = "red",
        bold = TRUE,
        font_size = 12
      ),
      kableExtra::cell_spec(x$`rhats > 1.05`, "html")
    )


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

# prints tables with fit coefficients and additional
html_format_coef_table <- function(x, fill, highlight){

  if (!is.null(x$Rhat)){
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
  }

  x_kbl <-
    kableExtra::kbl(
      x,
      row.names = TRUE,
      escape = FALSE,
      format = "html",
      digits = 3
    )

  if (highlight)
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

## helper modified from brms to get fit summary from draws object
draw_summary <- function(draws, variables, probs, robust) {

  .quantile <- function(x, ...) {
    qs <- posterior::quantile2(x, probs = probs, ...)
    prob <- probs[2] - probs[1]
    names(qs) <- paste0(c("l-", "u-"), prob * 100, "% CI")
    return(qs)
  }
  draws <- posterior::subset_draws(draws, variable = variables)
  measures <- list()
  if (robust) {
    measures$Estimate <- median
    # if (mc_se) {
    #   measures$MCSE <- posterior::mcse_median
    # }
    # measures$Est.Error <- mad
  } else {
    measures$Estimate <- mean
    # if (mc_se) {
    #   measures$MCSE <- posterior::mcse_mean
    # }
    # measures$Est.Error <- sd
  }

  measures$quantiles <- .quantile
  measures$Rhat <- posterior::rhat
  measures$Bulk_ESS <- posterior::ess_bulk
  measures$Tail_ESS <- posterior::ess_tail

  out <- do.call(posterior::summarize_draws, c(list(draws), measures))
  out <- as.data.frame(out)
  rownames(out) <- out$variable
  out$variable <- NULL
  return(out)
}
