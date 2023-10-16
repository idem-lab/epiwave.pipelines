
prepare_delay_input <- function(
        dates, jurisdictions,
        delay_data = NULL,
        constant_delay = NULL # should be value input if using
        ) {

    if (is.null(delay_data) & is.null(constant_delay)) {
        stop("Must specify either delay data or constant delay value")
    }

    if (!is.null(constant_delay)) {
        delay_input <- matrix(constant_delay,
                              nrow = length(dates),
                              ncol = length(jurisdictions),
                              dimnames = list(as.character(dates),
                                              jurisdictions))
        if (!is.null(delay_data)) {
            warning("Warning: delay data provided but overriden by constant delay")
        }
    }

    if (!is.null(delay_data) & is.null(constant_delay)) {
        delay_input <- estimate_delays_from_data(
            delay_data, jurisdictions)
    }

    delay_input
}



# check impact of on or multiple jurisdictions
estimate_delays_from_data <- function(linelist,
                                      jurisdictions) {

    # symptom_onset_date_col_name <- rlang::enquo(symptom_onset_date_col_name)
    # default state_specific delays, fall back on national

    # date_onset refers to symptom onset > change to date_symptom_onset?

    national_exclusions <- tibble::tibble(
        state = "VIC",
        start = as.Date("2020-06-14"),
        end   = as.Date("2020-12-01")
    )

    absolute_min_records <- 100

    # get delays for locally-acquired infections, truncated differently for
    # forward/backward delays
    detection_delay_data <- linelist |>
        dplyr::mutate(delay = as.numeric(
            date_confirmation - date_onset)) |>
        dplyr::filter(import_status == "local",
                      !is.na(date_onset),
                      delay <= 42,
                      delay >= -5)

    # WHAT DATE IS THIS SUPPOSED TO BE??
    # linelist_date <- local_linelist$date_linelist[1]
    # delay_data_from_onset <- detection_delay_data |>
    #     dplyr::filter(date_onset <= (linelist_date - 15))

    # THIS SHOULD BE WITH THE ABOVE OBJECT
    state <- detection_delay_data$state
    date <- detection_delay_data$date_onset
    delay <- detection_delay_data$delay

    # account for right-truncation when tabulating
    # which date to tabulate by
    date_from <- date
    date_to <- date + delay
    date_tabulation <- "date_from"

    delay_data <- tibble::tibble(
        state = state,
        date_from = date_from,
        date_to = date_to,
        delay = delay)

    all_dates <- seq(
        min(delay_data$date_from),
        max(delay_data$date_to),
        by = 1)

    all_states <- jurisdictions#unique(delay_data$state)

    date_state <- tidyr::expand_grid(
        date = all_dates,
        state = all_states)

    # get the half-window size (number of days on either side of the target)
    absolute_max_window <- as.numeric(diff(range(all_dates)))

    min_window <- 7
    max_window <- 56

    min_window <- ceiling((min_window - 1) / 2)
    max_window <- floor((max_window - 1) / 2)

    min_records <- 500

    # for each confirmation date, run the algorithm on each date
    statewide <- date_state |>
        dplyr::group_by(date, state) |>
        dplyr::mutate(window = get_window_size(
            date,
            state,
            delay_data = delay_data,
            date_tabulation = date_tabulation,
            n_min = min_records,
            window_min = min_window,
            window_max = max_window),
            count = count_in_window(
                date,
                state,
                delay_data = delay_data,
                date_tabulation = date_tabulation,
                window = window),
            state_ecdf = delay_ecdf(
                date,
                state,
                window = window,
                delay_data = delay_data,
                date_tabulation = date_tabulation))

    revert_to_national <- TRUE
    if (revert_to_national) {

        if(!is.null(national_exclusions)){
            # fill in exclusion periods
            national_exclusions <- national_exclusions |>
                dplyr::mutate(start = as.Date(start),
                              end = as.Date(end),
                              start = tidyr::replace_na(start, min(all_dates)),
                              end = tidyr::replace_na(end, max(all_dates)))

            # remove the specified data for estimating the
            #national background distribution
            for (i in seq_len(nrow(national_exclusions))) {
                delay_data <- delay_data |>
                    dplyr::filter(!(
                        state == national_exclusions$state[i] &
                            date_from >= national_exclusions$start[i] &
                            date_to <= national_exclusions$end[i]))}}

        nationwide <- date_state |>
            # arbitrarily pick one set of dates
            dplyr::filter(state == "NSW") |>
            dplyr::select(-state) |>
            dplyr::group_by(date) |>
            dplyr::mutate(window = get_window_size(
                date,
                all_states,
                delay_data = delay_data,
                date_tabulation = date_tabulation,
                n_min = min_records,
                window_min = min_window,
                window_max = absolute_max_window),
                national_ecdf = delay_ecdf(
                    date,
                    all_states,
                    window = window,
                    delay_data = delay_data,
                    date_tabulation = date_tabulation))

        # for statewide, replace any invalid ecdfs with the national one

        state_ecdfs <- statewide |>
            dplyr::right_join(nationwide |>
                                  dplyr::select(-window)) |>
            dplyr::mutate(
                use_national = count < absolute_min_records,
                weight = pmin(1, count / min_records),
                weight = ifelse(use_national, 0, weight),
                ecdf = mapply(
                    FUN = weight_ecdf,
                    state_ecdf,
                    national_ecdf,
                    weight,
                    SIMPLIFY = FALSE)) |>
            dplyr::select(date, state, ecdf, weight, use_national)

    } else {
        state_ecdfs <- statewide |>
            dplyr::mutate(use_national = count < absolute_min_records,
                          weight = pmin(1, count / min_records),
                          weight = ifelse(use_national, 0, weight),
                          ecdf = state_ecdf) |>
            dplyr::select(date, state, ecdf, weight, use_national)
    }


    # state_ecdfs is list of timevarying delays (either
    # functions or vector of probabilities)

    delay_dist <- state_ecdfs |>
        tidyr::pivot_wider(id_cols = date,
                           names_from = state,
                           values_from = ecdf)
    delay_dist_mat <- as.matrix(delay_dist[-1])
    row.names(delay_dist_mat) <- as.character(delay_dist$date)

    delay_dist_mat
}

count_in_window <- function(target_date, states,
                            delay_data, window,
                            date_tabulation) {
    dates <- delay_data |>
        dplyr::filter(state %in% states) |>
        dplyr::pull(!!date_tabulation)
    diff <- abs(dates - target_date)
    in_window <- diff <= window
    sum(in_window)
}

get_window_size <- function(
        target_date,
        states,
        delay_data,
        date_tabulation,
        n_min = 500,
        window_min = 7,
        window_max = 42) {

    dates <- delay_data |>
        dplyr::filter(state %in% states) |>
        dplyr::pull(!!date_tabulation)

    # find the smallest window that yields the required number of counts
    for (window in window_min:window_max) {

        diff <- abs(dates - target_date)
        in_window <- diff <= window

        if (sum(in_window) >= n_min) {
            break()
        }
    }
    window
}

delay_ecdf <- function(target_date, states, window,
                       delay_data, date_tabulation) {

    data <- delay_data |>
        dplyr::filter(state %in% states)
    dates <- dplyr::pull(data, !!date_tabulation)
    delays <- dplyr::pull(data, delay)

    diff <- abs(dates - target_date)
    in_window <- diff <= window

    valid_delays <- delays[in_window]

    if (length(valid_delays) > 0) {
        distribution <- ecdf(valid_delays)
    } else {
        distribution <- NULL
    }

    list(distribution)

}

# calculate a weighted average ecdf out of two
#(weight is the probability of the first)
weight_ecdf <- function(ecdf_1, ecdf_2, weight) {

    if (is.null(ecdf_1) | weight == 0) {
        return(ecdf_2)
    }
    if (is.null(ecdf_2) | weight == 1) {
        return(ecdf_1)
    }

    e1 <- environment(ecdf_1)
    e2 <- environment(ecdf_2)

    # reconcile the xs
    x_1 <- e1$x
    x_2 <- e2$x

    x <- sort(unique(c(x_1, x_2)))

    # get the two CDFs
    y_1 <- ecdf_1(x)
    y_2 <- ecdf_2(x)

    # get the two pdfs
    pdf_1 <- diff(c(0, y_1))
    pdf_2 <- diff(c(0, y_2))

    # get a weighted average of them
    pdf <- pdf_1 * weight + pdf_2 * (1 - weight)

    # convert back to a CDF
    y <- cumsum(pdf)

    # rebuild an ecdf object, the slow way
    method <- 2L
    yleft <- 0
    yright <- 1
    f <- e1$f
    n <- e1$nobs
    rval <- function (v) {
        stats:::.approxfun(x, y, v, method, yleft, yright, f)
    }
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- attr(ecdf_1, "call")
    rval

}

