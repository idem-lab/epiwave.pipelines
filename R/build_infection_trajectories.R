build_infection_trajectories <- function (param,
                                          infection_days,
                                          fit,
                                          infection_delay_distribution,
                                          jurisdictions = NULL) {

    infection_estimates <- build_estimates_array(
        param,
        infection_days,
        fit,
        jurisdictions)

    n_jurisdictions <- dim(infection_estimates)[3]

    # build convolution matrix
    out <- lapply(1:n_jurisdictions, function(x) {
        conv_mat <- lowerGPreff:::get_convolution_matrix(
            infection_delay_distribution[, x],
            nrow(infection_delay_distribution))
        inf_compl <- colSums(conv_mat)
        cbind(data.frame(jurisdiction = jurisdictions[x],
                         date = as.character(infection_days),
                         inf_compl = inf_compl),
              as.data.frame(infection_estimates[,,x]))
    })

    names(out) <- jurisdictions
    out
}

build_reff_estimates <- function (param,
                                  infection_days,
                                  fit,
                                  jurisdictions = NULL) {

    reff_estimates <- build_estimates_array(
        param,
        infection_days,
        fit,
        jurisdictions)

    n_jurisdictions <- dim(reff_estimates)[3]

    out <- lapply(1:n_jurisdictions, function(x) {

        df <- reff_estimates[,,x] |>
            as.data.frame()
        mat <- t(apply(df, 1, quantile,
                       c(0.05, 0.25, 0.5, 0.75, 0.95),
                       na.rm = TRUE))
        out <- cbind(data.frame(jurisdiction = jurisdictions[x],
                                date = as.character(infection_days)),
                     mat,
                     as.data.frame(reff_estimates[,,x]))
        out[-c(1:4),]

    })

    names(out) <- jurisdictions
    out
}
