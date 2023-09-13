fit_model <- function (combined_model_objects,
                       n_chains = 4,
                       max_convergence_tries = 2,
                       warmup = 500,
                       init_n_samples = 1000,
                       iterations_per_step = 1000) {

    m <- model(combined_model_objects$infections,
               combined_model_objects$expected_cases,
               combined_model_objects$gp_lengthscale)

    #get stable inits
    init <- generate_valid_inits(model = m,
                                 chains = n_chains,
                                 max_tries = 1000
    )

    # first pass at model fitting
    draws <- mcmc(
        m,
        sampler = hmc(Lmin = 25, Lmax = 30),
        chains = n_chains,
        warmup = warmup,
        n_samples = init_n_samples,
        initial_values = init,
        one_by_one = TRUE
    )

    # if it did not converge, try extending it a bunch more times
    finished <- converged(draws)
    tries <- 1
    while(!finished & tries < max_convergence_tries) {
        draws <- greta::extra_samples(
            draws,
            iterations_per_step,
            one_by_one = TRUE
        )
        tries <- tries + 1
        finished <- converged(draws)
    }

    # warn if we timed out before converging successfully
    if (tries == max_convergence_tries) {
        warning("sampling did not converge according to benchmarks")
    }

    return(draws)
}