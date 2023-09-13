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
