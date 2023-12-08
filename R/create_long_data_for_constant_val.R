create_long_data_for_constant_val <- function (
        dates,
        jurisdictions,
        constant_val,
        col_name) {

    long_unique <- expand.grid(date = dates,
                               jurisdiction = jurisdictions)

    long_combined <- cbind(
        long_unique,
        tibble::tibble(
            !!col_name := rep(constant_val, nrow = nrow(clean_pcr))
        )
    )
    long_combined
}
