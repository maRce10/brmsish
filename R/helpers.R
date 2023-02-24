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

# taken from brms:::escape_all
# .escape_all <- function (x)
# {
#   specials <- c(".", "*", "+", "?", "^", "$", "(", ")", "[",
#                 "]", "|")
#   for (s in specials) {
#     x <- gsub(s, paste0("\\", s), x, fixed = TRUE)
#   }
#   x
# }
#
# # taken from brms:::regex_or
# .regex_or <- function (x, escape = FALSE)
# {
#   if (escape) {
#     x <- .escape_all(x)
#   }
#   paste0("(", paste0("(", x, ")", collapse = "|"), ")")
# }
#
# # taken from brms:::valid_dpars
# .valid_dpars <- function (family, type = NULL, ...)
# {
#   if (!length(family)) {
#     if (is.null(type)) {
#       return("mu")
#     }
#     else {
#       return(NULL)
#     }
#   }
#   family <- .validate_family(family)
#   info <- paste0(.usc(type, "suffix"), "dpars")
#   .family_info(family, info, ...)
# }
#
# # taken from brms:::.family_info
# .family_info <- function (x, y, ...)
# {
#   x <- as_one_character(x)
#   y <- as_one_character(y)
#   if (y == "family") {
#     return(x)
#   }
#   if (!nzchar(x)) {
#     return(NULL)
#   }
#   info <- get(paste0(".family_", x))()
#   if (y == "link") {
#     out <- info$links[1]
#   }
#   else {
#     info$links <- NULL
#     out <- info[[y]]
#   }
#   out
# }
#
# # taken from brms:::validate_family
# .validate_family <- function (family, link = NULL, threshold = NULL)
# {
#   if (is.function(family)) {
#     family <- family()
#   }
#   if (!is(family, "brmsfamily")) {
#     if (is.family(family)) {
#       link <- family$link
#       family <- family$family
#     }
#     if (is.character(family)) {
#       if (is.null(link)) {
#         link <- family[2]
#       }
#       family <- .brmsfamily(family[1], link = link)
#     }
#     else {
#       stop2("Argument 'family' is invalid.")
#     }
#   }
#   if (is_ordinal(family) && !is.null(threshold)) {
#     threshold <- match.arg(threshold, c("flexible", "equidistant"))
#     family$threshold <- threshold
#   }
#   family
# }
#
# # taken from  brms:::usc
# .usc <- function (x, pos = c("prefix", "suffix"))
# {
#   pos <- match.arg(pos)
#   x <- as.character(x)
#   if (!length(x))
#     x <- ""
#   if (pos == "prefix") {
#     x <- ifelse(nzchar(x), paste0("_", x), "")
#   }
#   else {
#     x <- ifelse(nzchar(x), paste0(x, "_"), "")
#   }
#   x
# }
#
#
# .brmsfamily <- function (family, link = NULL, slink = link, threshold = "flexible",
#                          refcat = NULL, bhaz = NULL, ...)
# {
#   family <- tolower(as_one_character(family))
#   aux_links <- list(...)
#   pattern <- c("^normal$", "^zi_", "^hu_")
#   replacement <- c("gaussian", "zero_inflated_", "hurdle_")
#   family <- rename(family, pattern, replacement, fixed = FALSE)
#   ok_families <- lsp("brms", pattern = "^\\.family_")
#   ok_families <- sub("^\\.family_", "", ok_families)
#   if (!family %in% ok_families) {
#     stop2(family, " is not a supported family. Supported ",
#           "families are:\n", collapse_comma(ok_families))
#   }
#   family_info <- get(paste0(".family_", family))()
#   ok_links <- family_info$links
#   family_info$links <- NULL
#   if (!is.character(slink)) {
#     slink <- deparse(slink)
#   }
#   if (!slink %in% ok_links) {
#     if (is.character(link)) {
#       slink <- link
#     }
#     else if (!length(link) || identical(link, NA)) {
#       slink <- NA
#     }
#   }
#   if (length(slink) != 1L) {
#     stop2("Argument 'link' must be of length 1.")
#   }
#   if (is.na(slink)) {
#     slink <- ok_links[1]
#   }
#   if (!slink %in% ok_links) {
#     stop2("'", slink, "' is not a supported link ", "for family '",
#           family, "'.\nSupported links are: ", collapse_comma(ok_links))
#   }
#   out <- list(family = family, link = slink, linkfun = function(mu) link(mu,
#                                                                          link = slink), linkinv = function(eta) inv_link(eta,
#                                                                                                                          link = slink))
#   out[names(family_info)] <- family_info
#   class(out) <- c("brmsfamily", "family")
#   all_valid_dpars <- c(valid_dpars(out), valid_dpars(out, type = "multi"))
#   for (dp in all_valid_dpars) {
#     alink <- as.character(aux_links[[paste0("link_", dp)]])
#     if (length(alink)) {
#       alink <- as_one_character(alink)
#       valid_links <- links_dpars(dp)
#       if (!alink %in% valid_links) {
#         stop2("'", alink, "' is not a supported link ",
#               "for parameter '", dp, "'.\nSupported links are: ",
#               collapse_comma(valid_links))
#       }
#       out[[paste0("link_", dp)]] <- alink
#     }
#   }
#   if (is_ordinal(out$family)) {
#     thres_options <- c("flexible", "equidistant", "sum_to_zero")
#     out$threshold <- match.arg(threshold, thres_options)
#   }
#   if (conv_cats_dpars(out$family)) {
#     if (!has_joint_link(out$family)) {
#       out$refcat <- NA
#     }
#     else if (!is.null(refcat)) {
#       allow_na_ref <- !is_logistic_normal(out$family)
#       out$refcat <- as_one_character(refcat, allow_na = allow_na_ref)
#     }
#   }
#   if (is_cox(out$family)) {
#     if (!is.null(bhaz)) {
#       if (!is.list(bhaz)) {
#         stop2("'bhaz' should be a list.")
#       }
#       out$bhaz <- bhaz
#     }
#     else {
#       out$bhaz <- list()
#     }
#     if (is.null(out$bhaz$df)) {
#       out$bhaz$df <- 5L
#     }
#     if (is.null(out$bhaz$intercept)) {
#       out$bhaz$intercept <- TRUE
#     }
#   }
#   out
# }
#
# .as_one_character <- function (x, allow_na = FALSE)
# {
#   s <- substitute(x)
#   x <- as.character(x)
#   if (length(x) != 1L || anyNA(x) && !allow_na) {
#     s <- deparse_combine(s, max_char = 100L)
#     stop2("Cannot coerce '", s, "' to a single character value.")
#   }
#   x
# }
#
# .is.family <- function (x)
# {
#   inherits(x, "family")
# }
#
# .is_ordinal <- function (family)
# {
#   "ordinal" %in% .family_info(family, "specials")
# }
