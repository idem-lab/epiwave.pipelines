estimate_reff <- function (infections_timeseries,
                           generation_interval_mass_fxns) {

    convolution_matrix <- get_convolution_matrix(
        generation_interval_mass_fxns, nrow(infections_timeseries))
    infectiousness <- convolution_matrix %*% infections_timeseries

    reff <- infections_timeseries / (infectiousness + .Machine$double.eps)

    greta_arrays <- module(
        infectiousness,
        reff)

    return(greta_arrays)
}
