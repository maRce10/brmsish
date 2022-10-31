# check brmsfit models save as RDS files. Returns a data frame with parameters used to fit models. Can be used to make sure all models were run with the same parameters (e.g. before combining models)

check_rds_models <- function(path = ".", models = list.files(path = path, pattern = ".RDS$", ignore.case = TRUE, full.names = TRUE), parallel = 1, pb = TRUE){

    # run loop over models
    model_table_list <- pblapply_brmsish_int(X = models,cl = parallel, pbar = pb ,function(x){

        # set null objects to avoid errors when checking in CRAN
       Parameter <- NULL

        # read model
        model <- try(readRDS(x), silent = TRUE)

        if (!methods::is(model, 'try-error') & methods::is(model, 'brmsfit')){

        # get fit
        fit <- model$fit

        # get priors
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

        # get diagnostic quantities
        np <- brms::nuts_params(model)

        # put parameters on a table
        model_table <-
            data.frame(
                model = basename(x),
                b_prior,
                sd_prior,
                iterations = fit@stan_args[[1]]$iter,
                chains = length(attr(fit, "stan_args")),
                thinning = fit@stan_args[[1]]$thin,
                warmup = fit@stan_args[[1]]$warmup,
                n_parameters = length(names(fit@sim$samples[[1]])),
                diverg_transitions = sum(subset(np, Parameter == "divergent__")$Value),
                `rhats > 1.05` = sum(stats::na.omit(brms::rhat(model)) > 1.05),
                min_bulk_ESS = min( summary(model)$fixed$Bulk_ESS),
                min_tail_ESS = min( summary(model)$fixed$Tail_ESS),
                seed = fit@stan_args[[1]]$seed,
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

    if (length(model_table_list) == 0) cat("no readable models found:") else
    model_parameters <- do.call(rbind, model_table_list)

    return(model_parameters)

    }
