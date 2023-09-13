calculate_n_days_infection <- function(delay_list) {

    #infection burnin in days, defined by the max notification delay at the start
    max_burnin <- do.call(max, lapply(delay_list, function(df)
        max(which(df[[1]](0:28) != 0)))) # how long is longest possible delay
#df[[1]] pulls first step function (delay on first day)
    # outputs vector of probabilities of infection being delayed between 0 and 28 days

    max_days <- do.call(max, lapply(delay_list, function(df) length(df)))
# delays objects are list of prob mass functions?

    n_days_infection <- max_days + max_burnin
# do we want to also export the actual starting date?
    return(n_days_infection)
}

# a <- data.frame(c = 1:4, d = 1:4)
# b <- data.frame(c = 1:10, d = 1:10)

# max <- function (delay_list) {
#     do.call(max, lapply(delay_list, function(df) nrow(df)))
# }
# testfun(delay_list = list(a,b))

