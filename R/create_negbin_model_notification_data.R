create_negbin_model_notification_data <- function(
        gp_model_objects,
        timeseries_data,
        timevarying_delay_distribution,
        timevarying_proporion,
        n_days_infection,
        n_jurisdictions) {

    # fudge notification delay towards the end to match infection timepoints
    timevarying_delay_distribution[length(timevarying_delay_distribution):n_days_infection] <-
        timevarying_delay_distribution[length(timevarying_delay_distribution)]

    infections_timeseries <- gp_model_objects$infections_timeseries

    # convolve cases
    expected_cases <- convolve(timeseries_matrix = infections_timeseries,
                               mass_functions = timevarying_delay_distribution,
                               proportion = timevarying_proporion) # what does this look like when it's a matrix?

    # negative binomial likelihood for number of cases
    sqrt_inv_size <- normal(0, 0.5,
                            truncation = c(0, Inf),
                            dim = n_jurisdictions)
    sqrt_inv_size <- sweep(greta::zeros(n_days_infection,
                                        n_jurisdictions),
                           2, sqrt_inv_size,
                           FUN = "+")

    #calculate(sqrt_inv_size,nsim = 1)

    size <- 1 / sqrt(sqrt_inv_size)
    #calculate(size,nsim = 1)
    prob <- 1 / (1 + expected_cases / size)

    #observed days in the infection timeseries
    n_burnin <- n_days_infection - nrow(timeseries_data)
    # this line will also be specific to the shape of the input data
    obs_idx <- (n_burnin + 1) : (n_days_infection)

    expected_cases_obs <- expected_cases[obs_idx, ]
    size_obs <- size[obs_idx,]
    prob_obs <- prob[obs_idx,]
    #take out the extra infection days

    #define likelihood
    distribution(timeseries_data) <- greta::negative_binomial(
        size_obs,
        prob_obs)

    ## HOW to add new arrays to the gpmodobjects starting point.
    # do they need to be added here? or else just export the new ones?
    greta_arrays <- module(
        expected_cases,
        expected_cases_obs,
        size,
        prob,
        size_obs,
        prob_obs
    )
    # gp_model_objects[[length(gp_model_objects) + 1]] <- module(size)
    # gp_model_objects
    return(greta_arrays)
}
