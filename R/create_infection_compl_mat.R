create_infection_compl_mat <- function (convolution_matrices,
                                        jurisdictions) {

    with_named_jurisdictions <- setNames(convolution_matrices, #with_named_dates,
                                         jurisdictions)

    # we want proportion infections on that day that would be observed
    # we don't have
    # because this is centered on notification dates not infection dates
    # it tells us amount of cases notified to that date that we have a chance of observing
    # proportions of cases observable...
    inf_compl_prob_mat <- do.call(cbind, lapply(with_named_jurisdictions, colSums))
    inf_compl_prob_mat
}
