calculate_days_infection <- function(delay_list) {

    earliest_date <- as.Date(do.call(min, lapply(delay_list, function(df) min(rownames(df)))))

    #infection burnin in days, defined by the max notification delay at the start
    max_burnin <- do.call(max, lapply(delay_list, function(df)
        max(which(df[[1]](0:28) != 0)))) # how long is longest possible delay
    # outputs vector of probabilities of infection being delayed between 0 and 28 days

    latest_date <- as.Date(do.call(max, lapply(delay_list, function(df) max(rownames(df)))))

    infection_seq <- seq(
        from = earliest_date - max_burnin,
        to = latest_date,
        by = 'day')

    return(infection_seq)
}

# a <- data.frame(c = 1:4, d = 1:4)
# b <- data.frame(c = 1:10, d = 1:10)

# max <- function (delay_list) {
#     do.call(max, lapply(delay_list, function(df) nrow(df)))
# }
# testfun(delay_list = list(a,b))

