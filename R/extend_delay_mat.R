extend_delay_mat <- function (delay_matrix,
                              infection_days,
                              incubation_period) {

    delay_and_incubation <- apply(
        delay_matrix,
        c(1,2), construct_delays,
        ecdf2 = incubation_period_distribution,
        output = "probability",
        stepfun_output = TRUE)

    earliest_date <- as.Date(min(rownames(delay_and_incubation)))
    latest_date <- as.Date(max(rownames(delay_and_incubation)))

    extra_left <- as.character(earliest_date - min(PCR_infection_days))
    extra_right <- as.character(max(PCR_infection_days) - latest_date)

    delay_matrix_ext_left <- rbind(
        do.call("rbind",
                replicate(extra_left,
                          delay_and_incubation[1,],
                          simplify = FALSE)),
        delay_and_incubation)

    delay_matrix_ext_both <- rbind(
        delay_matrix_ext_left,
        do.call("rbind",
                replicate(extra_right,
                          delay_and_incubation[nrow(delay_and_incubation),],
                          simplify = FALSE)))

    rownames(delay_matrix_ext_both) <- as.character(infection_days)

    delay_matrix_ext_both
}

