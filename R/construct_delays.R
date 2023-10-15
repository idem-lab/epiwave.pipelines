#' Given a date range to compute delays over, and an empirical cdf, compute the probability of delay
#' having that length. Outputs either the probability or the cumulative density, as either a data
#' frame lined up against days-of-day or as a step function in the style of R native ecdf function.
#' Optionally, two cdfs can be used as input, in which case the probability and cdf of the additive
#' combined delay are calculated.
#'
#' @param delay_range
#' @param ecdf1
#' @param ecdf2
#' @param output
#' @param stefun_output
#'
#' @return
#' @export
#'
#' @examples
construct_delays <- function(ecdf1,
                             ecdf2 = NULL,
                             delay_range = c(-3, 28),
                             output = c("probability", "cumulative density"),
                             stefun_output = FALSE) {

if(is.list(ecdf1)) {
    ecdf1 <- ecdf1[[1]] }

    # days of delay
    days <- seq(delay_range[1], delay_range[2])

    if (is.null(ecdf2)) {
        # get discretised probability
        p <- ecdf1(days + 1) - ecdf1(days)
    } else {
        # get discretised probabilities for both cdfs
        p1 <- ecdf1(days + 1) - ecdf1(days)
        p2 <- ecdf2(days + 1) - ecdf2(days)

        p1 <- approxfun(days,p1, rule = 2)
        p2 <- approxfun(days,p2, rule = 2)

        #compute combined p
        p <- tidyr::expand_grid(
            x = days,
            z = days) |>
            dplyr::filter(x - z >= 0) |>
            dplyr::group_by(x) |>
            dplyr::summarise(p = sum(p1(z + 1) * p2(x - z + 1))) |>
            dplyr::pull(p)
    }

    #remove negative delay prob since notification cannot precede infection assuming that some
    #preventive measure has taken place once the "would-be" infectee is notified
    p[days < 0] <- 0
    #remove extremely long but unlikely delays, ususally a result of
    #parametric specification of the input cdf --- need to check if this is
    #right
    p[days > delay_range[2]] <- 0
    #normalise remaining probs
    p <- p / sum(p)

    if (output == "probability") {
        if (stefun_output) {
            return(approxfun(days,p,
                             rule = 2,
                             method = "constant",
                             yleft = 0,
                             yright = 0,
                             f = 0)
                   )
        } else {
            return(tibble::tibble(days,p))
        }
    } else {
        if (stefun_output) {
            #should we allow this? this literally just spits back out the input ecdf
            return(make_ecdf(p, days))
        } else {
            return(tibble::tibble(days, cumsum(p)))
        }
    }


}
