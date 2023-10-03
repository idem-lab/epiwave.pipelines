create_model_notification_data <- function(
        infections_timeseries,
        timevarying_delay_dist,
        timevarying_proportion,
        timeseries_data,
        model_likelihood = 'negative_binomial') {

    n_days_infection <- nrow(infections_timeseries)
    n_jurisdictions <- ncol(infections_timeseries)

    # fudge notification delay towards the end to match infection timepoints
    replication_times <- n_days_infection - nrow(timevarying_delay_dist)
    timevarying_delay_dist_ext <- rbind(
        timevarying_delay_dist,
        do.call("rbind",
                replicate(replication_times,
                          timevarying_delay_dist[nrow(timevarying_delay_dist),],
                          simplify = FALSE)))

    convolution_matrices <- lapply(1:n_jurisdictions, function(x)
        get_convolution_matrix(timevarying_delay_dist_ext[, x],
                               n_days_infection))

    expected_cases <- do.call(cbind, lapply(1:n_jurisdictions, function(x)
        convolution_matrices[[x]] %*% infections_timeseries[, x] * timevarying_proportion[, x]))

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
    if (model_likelihood == 'negative_binomial') {
        distribution(timeseries_data) <- greta::negative_binomial(
            size_obs,
            prob_obs)
    }

    greta_arrays <- module(
        convolution_matrices,
        expected_cases,
        expected_cases_obs,
        size,
        prob,
        size_obs,
        prob_obs
    )
    # infection_model_objects[[length(infection_model_objects) + 1]] <- module(size)
    # infection_model_objects
    return(greta_arrays)
}
