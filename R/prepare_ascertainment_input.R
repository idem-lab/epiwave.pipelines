# data usually pertains to behaviour around symptoms. make sure clear
# to user that CAR data should be related to notification not infection

prepare_ascertainment_input <- function(
        dates, jurisdictions,
        ascertainment_estimate = NULL,
        constant_ascertainment = NULL, # should be value input if using
        test_type = c("PCR", "RAT")) {

    # add option here likelihood for seroprevalence data to help constrain prior

    #sanity check for ascertainment specification
    if (is.null(ascertainment_estimate) &
        is.null(constant_ascertainment)) {

        stop("Provide either a constant ascertainment value or ascertainment
             estimate")
    }
    if (!is.null(constant_ascertainment)) {

        warning("Warning: using a constant assumed ascertainment value
                for every time point and jurisdiction!")

        ascertainment_input <- matrix(constant_ascertainment,
                                      nrow = length(dates),
                                      ncol = length(jurisdictions),
                                      dimnames = list(as.character(dates),
                                                      jurisdictions))
        if (!is.null(ascertainment_estimate)) {

            warning("Warning: ascertainment estimates provided but overriden by
                assumed constant ascertainment")
        }
    }
    if (!is.null(ascertainment_estimate) &
        is.null(constant_ascertainment)) {

        ascertainment_input <- get_CAR_from_surveys(
            dates, jurisdictions, ascertainment_estimate, test_type)
    }

    ascertainment_input

}


