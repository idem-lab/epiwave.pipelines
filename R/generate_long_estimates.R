#' Generate long estimates table
#'
#' @param param greta array for selected parameter
#' @param fit greta model draws
#' @param dates infection dates sequence
#' @param jurisdictions vector of jurisdictions
#'
#' @importFrom greta calculate
#' @importFrom dplyr mutate row_number
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_replace_all
#'
#' @return Long data frame with all draws and full date sequence
#' @export
generate_long_estimates <- function (param,
                                     fit,
                                     dates,
                                     jurisdictions) {

    sims <- greta::calculate(param,
                             values = fit,
                             nsim = 1000)[[1]]
    dates_df <- data.frame(ymd = dates,
                           date_id = as.character(1:length(dates)))

    df_list <- lapply(1:length(jurisdictions),
                      function (x) {
                          dat <- sims[ , -(1:5), x] |>
                              as.data.frame() |>
                              dplyr::mutate(.draw = dplyr::row_number(),
                                            jur = as.factor(jurisdictions[x]))})

    df <- do.call(rbind, df_list) |>
        tidyr::pivot_longer(-c(.draw, jur),
                            names_to = 'date_id',
                            values_to = 'val') |>
        dplyr::mutate(date_id = as.numeric(
            stringr::str_replace_all(date_id, 'V', ''))) |>
        merge(dates_df)

    df
}

