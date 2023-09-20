create_timevarying_delay_distribution <- function(
        notification_matrix, # in real one shouldn't have this as input
        incubation_period_distribution) {

    dates <- seq_len(nrow(notification_matrix))
    # # set date names
    # if (is.null(dates)) {
    #     dates <- rownames(notification_matrix)
    # } else if (length(dates) != nrow(notification_matrix)) {
    #     stop("Error: supplied date labels has different length from number of rows in notification matrix!")
    # } else {
    #     dates <- dates
    # }

    n_days <- nrow(notification_matrix) # why not also length(dates) ?

    #fake some delay data
    set.seed(2023-06-26)

    notification_delay_data <- tibble::tibble(
        days = sample(seq_along(dates), length(dates)*100, replace = TRUE),
        delay = ceiling(rnorm(length(dates)*100, mean = 4, sd = 0.5)) #rep(3,n_days*100)
    )

    # fake delay ecdf
    notification_delay_function_raw <- estimate_delays(
        notification_delay_data = notification_delay_data,
        time_varying = TRUE #currently testing if time-fixed delay help with convergence
    )

    #quick check delay behaviour
    # lapply(1:5,FUN = function(x){ notification_delay_function_raw$delay_ecdf[x][[1]](1:5)})
    # add functions to combine an incubation period distribution ( from
    # literature) with an onset to notification distribution ( from linelist
    # data)

    # combine the incubation and notification delays
    timevarying_delay_distribution <- lapply(
        notification_delay_function_raw$delay_ecdf,
        FUN = function(x) {
            delay_constructor(ecdf1 = x,
                              ecdf2 = incubation_period_distribution,
                              output = "probability",
                              stefun_output = TRUE)
        }
    )

    return(timevarying_delay_distribution)
}

#' Takes a data frame of "notification_delay_data", and returns a delay ecdf
#' across all time or a list of ecdfs for each time step. The
#' notification_delay_data should have a column of days (or other time
#' intervals) and a column of FORWARD delays from the corresponding days. This
#' is a fudge function to quickly put together delay data in the format usable
#' for the rest of the code, the real delay data will have more missingness each
#' day, and therefore need to be estimated with a rolling time window, in a more
#' complicated way
#'
#' @param notification_delay_data
#' @param time_varying
#'
#' @return
#' @export
#'
#' @examples
estimate_delays <- function(notification_delay_data,
                            time_varying = TRUE) {

    if (time_varying) {
        delays <- notification_delay_data |>
            dplyr::group_by(days) |>
            dplyr::summarise(delay_ecdf = list(ecdf(delay)))
    } else {
        delays <- notification_delay_data |>
            dplyr::group_by(days) |>
            dplyr::summarise(delay_ecdf = list(ecdf(notification_delay_data$delay)))
    }

    # #save choice for time varying or not
    # return(list(
    #     delays = delays,
    #     time_varying = time_varying
    #     )
    # )
    #save just the delays
    return(delays)
}


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
delay_constructor <- function(delay_range = c(-3, 28),
                              ecdf1,
                              ecdf2 = NULL,
                              output = c("probability", "cumulative density"),
                              stefun_output = FALSE) {


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

        #remove negative delay prob since notification cannot precede infection assuming that some
        #preventive measure has taken place once the "would-be" infectee is notified
        p[days < 0] <- 0
        p <- p / sum(p)

    }

    if (output == "probability") {
        if (stefun_output) {
            return(approxfun(days,p, rule = 2))
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
