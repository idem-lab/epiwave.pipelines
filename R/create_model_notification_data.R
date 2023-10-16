create_model_notification_data <- function(
        infections_timeseries,
        full_infection_dates,
        observed_infection_dates,
        timevarying_delay_dist_ext,
        timevarying_proportion,
        timeseries_data,
        model_likelihood = 'negative_binomial') {

    n_days_infection <- nrow(infections_timeseries)
    n_jurisdictions <- ncol(infections_timeseries)

    #get indices for subset of infection days that end up in observed/forecast data
    obs_infection_idx <- which(full_infection_dates %in% observed_infection_dates)
    #subset infection timeseries to these indices
    infection_obs <- infections_timeseries[obs_infection_idx,]
    n_days_infection_obs <- length(obs_infection_idx)

    convolution_matrices <- lapply(1:n_jurisdictions, function(x)
        get_convolution_matrix(timevarying_delay_dist_ext[, x],
                               n_days_infection_obs))
    #compute expected cases of the same length
    #note not all of these dates would have been observed
    expected_cases <- do.call(cbind,
                              lapply(
        1:n_jurisdictions, function(x)
        {convolution_matrices[[x]] %*% infection_obs[, x] * timevarying_proportion[, x]}
                                    )
                            )

    # negative binomial likelihood for number of cases
    sqrt_inv_size <- normal(0, 0.5,
                            truncation = c(0, Inf),
                            dim = n_jurisdictions)
    sqrt_inv_size <- sweep(greta::zeros(n_days_infection_obs,
                                        n_jurisdictions),
                           2, sqrt_inv_size,
                           FUN = "+")

    size <- 1 / sqrt(sqrt_inv_size)
    prob <- 1 / (1 + expected_cases / size)

    #observed days in the data itself
    obs_data_idx <- which(observed_infection_dates %in% as.Date(rownames(timeseries_data)))
    #basically this backs out burn in and burn out days
    size_obs <- size[obs_data_idx,]
    prob_obs <- prob[obs_data_idx,]
    #this subsets to actually observed data dates
    #expected_cases_obs <- expected_cases[obs_data_idx,]
    infection_match_data <- infection_obs[obs_data_idx,]

    #define likelihood
    if (model_likelihood == 'negative_binomial') {

        #realistically you only care about forecast if using case data, thus
        #this bit is conditional on neg binom likelihood
        #get idx for forecast days
        forecast_data_idx <- which(observed_infection_dates > max(as.Date(rownames(timeseries_data)))
        )
        #proposal for a forecast bit
        size_forecast <- size[forecast_data_idx,]
        prob_forecast <- prob[forecast_data_idx,]
        forecast_cases <- greta::negative_binomial(
            size_forecast,
            prob_forecast)
        #array-ify the data
        timeseries_data_array <- greta::as_data(timeseries_data)
        #stick in the forecast bit
        timeseries_data_array_with_forecast <- rbind(timeseries_data_array,
                                                     forecast_cases)

        #I guess you could make a version of infection time series with forecast
        #bit too but can't imagine that being useful at all
        distribution(timeseries_data_array) <- greta::negative_binomial(
            size_obs,
            prob_obs)
    }

    greta_arrays <- module(
        infection_obs,
        infection_match_data,
        size,
        prob,
        timeseries_data_array,
        timeseries_data_array_with_forecast
    )
    # infection_model_objects[[length(infection_model_objects) + 1]] <- module(size)
    # infection_model_objects
    return(greta_arrays)
}
