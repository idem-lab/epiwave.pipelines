build_estimates_array <- function (param,
                                   infection_days,
                                   fit,
                                   jurisdictions = NULL,
                                   nsim = 1000) {

    # calculate infection_timeseries from model
    sims <- greta::calculate(param,
                             values = fit,
                             nsim = nsim)[[1]]

    if (is.null(jurisdictions) & !is(sims, 'matrix')) {
        stop('List jurisdiction names to match modelled arrays')
    }
    if (is.null(jurisdictions) & is(sims, 'matrix')) {
        n_jurisdictions <- 1
        sims_array <- array(sims,
                            dim = c(nrow(sims),
                                    ncol(sims), 1))
    }
    if (!is.null(jurisdictions)) {
        n_jurisdictions <- length(jurisdictions)
        sims_array <- sims
    }

    transposed <- aperm(sims_array,
                        c(2, 1, 3))
    draws <- paste0('draw', 1:nsim)
    dimnames(transposed) <- list(
        NULL,
        draws,
        jurisdictions)

    transposed
}
