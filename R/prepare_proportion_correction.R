prepare_proportion_correction <- function(
        dates, # test type specific date sequence
        observable_infection_dates,
        timevarying_ascertainment,
        case_type_proportion = NULL,
        dow_correction = NULL
        # dates, jurisdictions,
        # proportion_estimate = NULL,
        # constant_ascertainment = NULL, # should be value input if using
        # test_type = c("PCR", "RAT")
) {

    # observed days in the obs data itself
    # basically this backs out extra left and right days
    obs_data_idx <- which(observable_infection_dates %in% dates)

    # update proportion arg with type proportion
    if (!is.null(case_type_proportion)) {
        case_type_proportion <- pmax(case_type_proportion, 1e-3)
        timevarying_ascertainment[obs_data_idx,] <-
            timevarying_ascertainment[obs_data_idx,] * case_type_proportion
    }

    # update proportion arg with dweek correction
    if (is.null(dow_correction)) {
        timevarying_proportion <- timevarying_ascertainment
    }
    if (!is.null(dow_correction)) {
        timevarying_proportion <- timevarying_ascertainment * dow_correction
    }


    out <- list(timevarying_proportion = timevarying_proportion)
    out

#
#
#
#     ##### these copied from CAR fxn
#     # add option here likelihood for seroprevalence data to help constrain prior
#
#     #sanity check for ascertainment specification
#     if (is.null(ascertainment_estimate) &
#         is.null(constant_ascertainment)) {
#
#         stop("Provide either a constant ascertainment value or ascertainment
#              estimate")
#     }
#     if (!is.null(constant_ascertainment)) {
#
#         warning("Warning: using a constant assumed ascertainment value
#                 for every time point and jurisdiction!")
#
#         ascertainment_input <- matrix(constant_ascertainment,
#                                       nrow = length(dates),
#                                       ncol = length(jurisdictions),
#                                       dimnames = list(as.character(dates),
#                                                       jurisdictions))
#         if (!is.null(ascertainment_estimate)) {
#
#             warning("Warning: ascertainment estimates provided but overriden by
#                 assumed constant ascertainment")
#         }
#     }
#     if (!is.null(ascertainment_estimate) &
#         is.null(constant_ascertainment)) {
#
#         ascertainment_input <- get_CAR_from_surveys(
#             dates, jurisdictions, ascertainment_estimate, test_type)
#     }
#
#     ascertainment_input

}
