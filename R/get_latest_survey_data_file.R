get_latest_survey_data_file <- function() {

    total_estimate <- list.files(
        'outputs',
        pattern = 'at_least_one_sym_states_central_smoothed_20',
        full.names = TRUE
    )

    file_dates <- total_estimate |>
        substr(50, 59) |>
        as.Date(format = '%Y-%m-%d')

    latest <- which.max(file_dates)
    latest_file <- total_estimate[latest]

    latest_file

}
