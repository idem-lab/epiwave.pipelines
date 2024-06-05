#' Expand constant value into long tibble
#'
#' @description The epiwave model functions expect data in a long format,
#'  which is structed to have a value for every unique date and jurisdiction
#'  pair. This function create a tibble of this structure out of a single
#'  value that should be replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param jurisdictions jurisdiction names
#' @param constant_val value to be replicated in each cell
#' @param col_name column name to contain the `constant value`
#'
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#'
#' @return long tibble with constant value replicated
#' @export
expand_constant_value <- function (dates,
                                   jurisdictions,
                                   constant_val,
                                   col_name) {

    long_unique <- expand.grid(date = dates,
                               jurisdiction = jurisdictions)
    long_combined <- long_unique |>
        dplyr::mutate(!!col_name := constant_val) |>
        tibble::tibble()

    long_combined
}
