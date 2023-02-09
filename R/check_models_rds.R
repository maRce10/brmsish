# check brmsfit models save as RDS files. Returns a data frame with parameters used to fit models. Can be used to make sure all models were run with the same parameters (e.g. before combining models)

check_rds_models <- function(path = ".", models = list.files(path = path, pattern = ".RDS$", ignore.case = TRUE, full.names = TRUE), parallel = 1, pb = TRUE, robust = TRUE, html = FALSE){

    # run loop over models
    model_table_list <- pblapply_brmsish_int(X = models, cl = parallel, pbar = pb ,function(x){

        # set null objects to avoid errors when checking in CRAN
       Parameter <- NULL

        # read model
        model <- try(readRDS(x), silent = TRUE)

        if (!methods::is(model, 'try-error') & methods::is(model, 'brmsfit')){

        # get diagnostic quantities
        iterations <- model$fit@sim$iter
        chains <- posterior::nchains(model)
        warmup <- model$fit@sim$warmup
        mod_formula <- as.character(model$formula[1])
        diverg_transitions <- sum(brms::nuts_params(model, pars = "divergent__")$Value)
        priors <- paste(apply(brms::prior_summary(model, all = FALSE)[, 2:1], 1, paste, collapse = "-"), collapse = "\n")
        seed <- model$fit@stan_args[[1]]$seed
        thinning <- model$fit@stan_args[[1]]$thin
        n_parameters <- length(names(model$fit@sim$samples[[1]]))

        # get draws
        vars <- posterior::variables(model)
        model <- posterior::as_draws_array(model, variable = vars)

        # get summary
        coef_table <- draw_summary(model, variables = vars, probs = c(0.025, 0.975), robust = robust)

        # put parameters on a table
        model_table <-
          data.frame(
            model = basename(x),
            priors = priors,
            formula = mod_formula,
            iterations = iterations,
            chains = chains,
            thinning = thinning,
            warmup = warmup,
            n_parameters = n_parameters,
            diverg_transitions = diverg_transitions,
            `rhats > 1.05` = sum(coef_table$Rhat > 1.05),
            min_bulk_ESS = min(coef_table$Bulk_ESS),
            min_tail_ESS = min(coef_table$Tail_ESS),
            seed = seed,
            check.names = FALSE
          )

        } else
            model_table <- NA

            return(model_table)

    })

    # which were not read
    which_data_frame <- sapply(model_table_list, is.data.frame)

    not_read <-basename(models)[!which_data_frame]

    if (length(not_read) > 0){
        cat("The following files cannot be read or are not 'brmsfit' models:")
    print(not_read)

    model_table_list <- model_table_list[which_data_frame]
     }

    if (length(model_table_list) == 0) cat("no readable models found") else
    model_parameters <- do.call(rbind, model_table_list)

    if (!html)
    return(model_parameters) else {

      html_tab <- html_format_model_table(model_parameters)

      return(html_tab)
    }

    }
