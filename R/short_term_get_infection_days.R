short_term_get_infection_days <- function (case_dat,
                                           delay_dat) {

    dates <- unique(case_dat$date)
    jurisdictions <- unique(case_dat$jurisdiction)
    delay_mat <- prepare_delay_input(
        dates, jurisdictions,
        constant_delay = list(delay_dat))
    calculate_days_infection(delay_mat)

}
calculate_days_infection <- function (delay_matrix) {

    # extra_left is defined by the max notification delay at the start
    extra_left <- max(which(delay_matrix[[1]](0:28) != 0))

    # extra_right defined by the max notification delay at the end
    extra_right <- max(which(delay_matrix[[nrow(delay_matrix)]](0:28) != 0))

    earliest_date <- as.Date(min(rownames(delay_matrix)))

    # outputs vector of probabilities of infection being delayed between 0 and 28 days
    latest_date <- as.Date(max(rownames(delay_matrix)))

    infection_seq <- seq(
        from = earliest_date - extra_left,
        to = latest_date + extra_right,
        by = 'day')

    return(infection_seq)
}

