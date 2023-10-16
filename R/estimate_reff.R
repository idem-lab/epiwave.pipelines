estimate_reff <- function (infections_timeseries,
                           generation_interval_mass_fxns) {

    # infectiousness <- convolve(timeseries_matrix = infections_timeseries,
    #                            mass_functions = generation_interval_mass_fxns,
    #                            proportion = 1)

    convolution_matrix <- get_convolution_matrix(
        generation_interval_mass_fxns, nrow(infections_timeseries))
    infectiousness <- convolution_matrix %*% infections_timeseries

    reff <- infections_timeseries / (infectiousness + .Machine$double.eps)

    greta_arrays <- module(
        infectiousness,
        reff)

    return(greta_arrays)
}
