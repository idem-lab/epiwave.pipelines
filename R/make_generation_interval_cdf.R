# the serial interval of Nishiura et al.
make_generation_interval_cdf <- function(
        generation_interval_distribution_data) {

    # generation interval distribution; use SI distribution from Nishiura et al.
    meanlog <- mean(generation_interval_distribution_data$param1)
    sdlog <- mean(generation_interval_distribution_data$param2)

    # this has to be kept inside the parent function so the meanlog and
    # sdlog pass correctly
    gi_cdf <- function(days) {
        plnorm(days, meanlog, sdlog)
    }

    gi_cdf

}

