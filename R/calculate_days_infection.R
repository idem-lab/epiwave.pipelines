calculate_days_infection <- function(delay_matrix) {

    earliest_date <- as.Date(
        min(
            rownames(delay_matrix)
            )
        )

    #infection burnin in days, defined by the max notification delay at the start
    max_burnin <- max(which(delay_matrix[[1]](0:28) != 0)
                      ) # how long is longest possible delay
    # outputs vector of probabilities of infection being delayed between 0 and 28 days
    latest_date <- as.Date(
        max(
            rownames(delay_matrix)
        )
    )
    #infection burnin out days, defined by the max notification delay at the end
    max_burnout <- max(which(delay_matrix[[nrow(delay_matrix)]](0:28) != 0)
    ) # how long is longest possible delay

    infection_seq <- seq(
        from = earliest_date - max_burnin,
        to = latest_date + max_burnout,
        by = 'day')

    return(infection_seq)
}


