create_infection_timeseries <- function(n_jurisdictions,
                                        n_days_infection,
                                        effect_type = "growth_rate") {

    # kernel hyperparams
    gp_lengthscale <- greta::lognormal(0, 1) #inverse_gamma(187/9,1157/18)
    gp_variance <- greta::normal(0, 4, truncation = c(0, Inf))
    gp_kernel <- greta.gp::mat52(gp_lengthscale, gp_variance)

    #define gp
    gp <- greta.gp::gp(
        x = seq_len(n_days_infection),
        kernel = gp_kernel,
        n = n_jurisdictions,
        tol = 1e-4
    )

    #compute infections from gp
    infections_timeseries <- compute_infections(
        log_effect = gp,
        effect_type = effect_type
    )

    greta_arrays <- module(
        gp,
        infections_timeseries,
        gp_lengthscale,
        gp_variance,
        gp_kernel
    )

    return(greta_arrays)
}
