summarise_linelist <- function (linelist,
                                import_status_option) {

    n_cases_line_list <- linelist |>
        dplyr::filter(import_status == import_status_option) |>
        dplyr::group_by(date_confirmation, date_onset, state, test_type) |>
        dplyr::summarise(n_cases = dplyr::n()) |>
        tidyr::pivot_wider(id_cols = c(date_confirmation, date_onset, state),
                           names_from = test_type,
                           values_from = n_cases,
                           values_fill = 0) |>
        dplyr::mutate(total = sum(PCR, RAT, na.rm = TRUE))

    n_cases_line_list

}
