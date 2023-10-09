
prepare_ascertainment_input <- function(
        dates, jurisdictions,
        ascertainment_estimate = NULL,
        assume_constant_ascertainment = NULL) { # should be value

    # add option here likelihood for seroprevalence data to help constrain prior

    #sanity check for ascertainment specification
    if (is.null(ascertainment_estimate) & is.null(assume_constant_ascertainment)) {
        warning("Warning: not assuming constant ascertainment but no ascertainment estimates provided, ascertainment rate with be fitted with a naive prior")
    } else if (!is.null(ascertainment_estimate) & !is.null(assume_constant_ascertainment)) {
        warning("Warning: ascertainment estimates provided but overriden by assumed constant ascertainment")
    } else if (is.null(ascertainment_estimate) & !is.null(assume_constant_ascertainment)) {
        ascertainment_input <- matrix(assume_constant_ascertainment,
                                      nrow = length(dates),
                                      ncol = length(jurisdictions),
                                      dimnames = list(dates, jurisdictions))
    } else if (!is.null(ascertainment_estimate) & is.null(assume_constant_ascertainment)) {
        ascertainment_input <- get_CAR_from_surveys(dates, jurisdictions, ascertainment_estimate)
    }

    ascertainment_input

}


# date_state_mat <- create_date_state_matrix()
# latest_survey_data <- get_latest_survey_data_file()
#
get_CAR_from_surveys <- function (dates, jurisdictions, latest_survey_data) {

    date_state <- dplyr::tibble(
        date = rep(dates, length(jurisdictions)),
        state = rep(jurisdictions, each = length(dates)))

    CAR_time_point <- latest_survey_data |>
        readr::read_csv() |>
        dplyr::select(state, date, fitted_point4) |>
        dplyr::rename('test_prob_given_symptom' = fitted_point4) |>
        dplyr::mutate(test_prob_given_infection = test_prob_given_symptom * 0.75 * 0.01)

    #join in the estimates
    CAR_smooth <- dplyr::left_join(tibble::tibble(date_state),
                                   CAR_time_point)

    #remove conditional on symptom column since it's linearly scaled from infection
    CAR_smooth <- CAR_smooth |>
        dplyr::select(-test_prob_given_symptom)

    last_perfection_ascertainment_date <- lubridate::as_date('2021_12_07')

    #hack in 0.75 for Delta period
    CAR_smooth <- CAR_smooth |>
        dplyr::mutate(test_prob_given_infection =
                          dplyr::case_when(
                              date <= last_perfection_ascertainment_date ~ 0.75,
                              TRUE ~ test_prob_given_infection))

    #Add quick drop off in Dec 2021
    CAR_smooth <- CAR_smooth |>
        dplyr::mutate(test_prob_given_infection =
                          dplyr::case_when(
                              date == as.Date('2021-12-14') ~ 0.33,
                              TRUE ~ test_prob_given_infection))

    #interpolate NAs and repeat forward last survey result
    CAR_smooth <- CAR_smooth |>
        dplyr::group_by(state) |>
        dplyr::mutate(test_prob_given_infection =
                          zoo::na.approx(test_prob_given_infection,
                                         na.rm = FALSE),
                      test_prob_given_infection =
                          zoo::na.locf(test_prob_given_infection,
                                       na.rm = FALSE,
                                       fromLast = FALSE)) |>
        dplyr::rename(CAR = test_prob_given_infection)

    #get matrix form
    CAR_smooth_mat <- CAR_smooth |>
        dplyr::select(date,state,CAR) |>
        tidyr::pivot_wider(values_from = CAR,
                           names_from = state) |>
        dplyr::select(-date) |>
        as.matrix()

    CAR_smooth_mat

}

