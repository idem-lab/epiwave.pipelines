
prepare_ascertainment_input <- function(
        assume_constant_ascertainment = FALSE,
        ascertainment_estimate = NULL,
        constant_ascertainment_fraction = 1) {

    # add option here likelihood for seroprevalence data to help constrain prior

    #sanity check for ascertainment specification
    if (is.null(ascertainment_estimate) & !assume_constant_ascertainment) {
        warning("Warning: not assuming constant ascertainment but no ascertainment estimates provided, ascertainment rate with be fitted with a naive prior")
    } else if (!is.null(ascertainment_estimate) & assume_constant_ascertainment) {
        warning("Warning: ascertainment estimates provided but overriden by assumed constant ascertainment")
    } else if (is.null(ascertainment_estimate) & assume_constant_ascertainment) {
        ascertainment_input <- constant_ascertainment_fraction
    } else if (!is.null(ascertainment_estimate) & !assume_constant_ascertainment) {
        # do something here to use the ascertainment estimates to inform prior
        #dummy function below
        ascertainment_input <- get_ascertainment_prior(ascertainment_estimate)
        # WHAT IS THE CONTENT OF THIS STEP?
    }

}
