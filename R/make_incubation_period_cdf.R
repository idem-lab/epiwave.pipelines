# takes estimates of incubation period from literature and make it a cdf
#' Note that in using incubation period for modelling we should not treat it as time-varying (due to
#' changes in the dominant variant). Maybe through vignette or some other documentation we should
#' make it clear that this is not a time-varying biological quantity and the change in dominant
#' variant should be modelled separately)
#'
#'
#' @param strain
#'
#' @return
#' @export
#'
#' @examples
make_incubation_period_cdf <- function(
        strain = c("WT",
                   "Alpha",
                   "Beta/Gamma",
                   "Delta",
                   "Omicron")) {

    strain <- match.arg(strain)

    if (strain == "Omicron") {

        days <- 0:28
        cum_density <- pweibull(days,shape = 1.5,scale = 3.6)
        # replaced these with a better one now
        # prob <- c(0.05, 0.5, 0.95)
        # estimate <- c(2.88, 3.42, 3.96)
    }

    cdf <- approxfun(days,cum_density)
    #cdf <- make_ecdf(prob, estimate)
    return(cdf)
}
