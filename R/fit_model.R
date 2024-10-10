#' Fit model
#'
#' @param model greta model object
#' @param n_chains number of chains
#' @param max_convergence_tries maximum number of times extra samples should be
#'  drawn if convergence isn't reached
#' @param warmup number of warmup samples
#' @param n_samples number of additional samples after warmup
#' @param n_extra_samples number of samples to be taken with `extra_samples()`
#' @param one_by_one whether to run TensorFlow MCMC code one iteration at a time,
#'  so that greta can handle numerical errors as 'bad' proposals.
#'
#' @importFrom greta extra_samples hmc mcmc
#'
#' @return fitted draws
#' @export
fit_model <- function (model,
                       n_chains = 4,
                       max_convergence_tries = 5,
                       warmup = 1000,
                       n_samples = 1000,
                       n_extra_samples = 1000,
                       one_by_one = TRUE) {

    # # get stable inits
    # init <- generate_valid_inits(model = model,
    #                              chains = n_chains,
    #                              max_tries = 1000)

    # first pass at model fitting
    draws <- greta::mcmc(
        model,
        sampler = greta::hmc(Lmin = 25, Lmax = 30),
        chains = n_chains,
        warmup = warmup,
        n_samples = n_samples,
        #initial_values = init,
        one_by_one = one_by_one
    )

    # if it did not converge, try extending it a bunch more times
    finished <- check_convergence(draws)
    tries <- 1
    while(!finished & tries < max_convergence_tries) {
        draws <- greta::extra_samples(
            draws,
            n_extra_samples,
            one_by_one = TRUE
        )
        tries <- tries + 1
        finished <- check_convergence(draws)
    }

    # warn if we timed out before converging successfully
    if (tries == max_convergence_tries) {
        warning("sampling did not converge according to benchmarks")
    }

    return(draws)
}

#' Check convergence
#'
#' @description Check whether the sampler has converged to our standards.
#'
#' @param draws fitted draws
#' @param max_r_hat maximum acceptable r hat
#' @param min_n_eff minimum acceptable n effective
#'
#' @return TRUE or FALSE on whether convergence was achieved per parameters
check_convergence <- function (draws,
                               max_r_hat = 1.1,
                               min_n_eff = 1000) {

    stats <- get_convergence_stats(draws)
    all(stats$r_hats < max_r_hat) &
        all(stats$n_eff >= min_n_eff)

}

#' Get convergence stats
#'
#' @param draws fitted draws
#'
#' @importFrom coda effectiveSize gelman.diag
#'
#' @return invisibly return list of max r hat and min n effective
get_convergence_stats <- function (draws) {

    r_hats <- coda::gelman.diag(draws, autoburnin = FALSE,
                                multivariate = FALSE)$psrf[, 1]
    n_eff <- coda::effectiveSize(draws)

    # sometimes n_eff underflows to 0 simply because the values being traced are
    # very small, so remove these (exactly 0 is not possible)
    n_eff <- n_eff[n_eff != 0]

    cat(sprintf("maximum R-hat: %.2f\nminimum n effective: %.2f",
                max(r_hats, na.rm = TRUE),
                min(n_eff)))

    result <- list(r_hats = r_hats,
                   n_eff = n_eff)

    invisible(result)
}

