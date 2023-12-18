#' Build table of parameter trajectory estimates
#'
#' @description Build an output table with all the simulations in order
#'  to pass along to the decision support engine.
#'
#' @param param model parameter of interest
#' @param infection_days full date sequence of infection timeseries
#' @param fit model fit
#' @param nsim number of simulations from model - must be <= samples taken
#' @param jurisdictions optional vector of jurisdictions
#'
#' @importFrom greta calculate
#' @importFrom methods is
#'
#' @return table with jurisdiction, date, and all draws for given parameter
#' @export
build_trajectories <- function (param,
                                infection_days,
                                fit,
                                nsim = 1000,
                                jurisdictions = NULL) {

    # calculate infection_timeseries from model
    sims <- greta::calculate(param,
                             values = fit,
                             nsim = nsim)[[1]]

    if (is.null(jurisdictions) & !methods::is(sims, 'matrix')) {
        stop('List jurisdiction names to match modelled arrays')
    }
    if (is.null(jurisdictions) & methods::is(sims, 'matrix')) {
        n_jurisdictions <- 1
        sims_array <- array(sims,
                            dim = c(nrow(sims),
                                    ncol(sims), 1))
    }
    if (!is.null(jurisdictions)) {
        n_jurisdictions <- length(jurisdictions)
        sims_array <- sims
    }

    estimates <- aperm(sims_array,
                        c(2, 1, 3))
    draws_id <- paste0('draw', 1:nsim)
    dimnames(estimates) <- list(
        NULL,
        draws_id,
        jurisdictions)

    n_jurisdictions <- dim(estimates)[3]

    # build output table
    out <- lapply(1:n_jurisdictions, function(x) {
        cbind(data.frame(jurisdiction = jurisdictions[x],
                         date = as.character(infection_days)),
              as.data.frame(estimates[,,x]))
    })

    names(out) <- jurisdictions
    out
}
