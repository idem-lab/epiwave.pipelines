create_dow_correction_objects <- function (
        observable_infection_dates,
        n_jurisdictions,
        dataID) {

    # get day of week index
    # favour base approach over lubridate for manual declaration of factor level
    # so that we know which index is which day of the week
    dweek <- weekdays(observable_infection_dates)
    dweek <- as.integer(
        factor(
            dweek,
            levels = c("Monday",
                       "Tuesday",
                       "Wednesday",
                       "Thursday",
                       "Friday",
                       "Saturday",
                       "Sunday")
        )
    )

    # prior for dweek correction
    dow_alpha <- greta::normal(1, 1, truncation = c(0, Inf), dim = c(1, 7))
    dow_dist <- greta::dirichlet(dow_alpha, n_realisations = n_jurisdictions)

    # normalise multiplier to average to 1
    dow_weights <- dow_dist * 7

    # match weight to date by state matrix
    dow_correction <- t(dow_weights[, dweek])

    dow_greta_arrays <- list(
        dow_weights,
        dow_correction)

    names(dow_greta_arrays) <- c(
        paste0(dataID, '_dow_weights'),
        paste0(dataID, '_dow_correction'))

    return(dow_greta_arrays)
}
