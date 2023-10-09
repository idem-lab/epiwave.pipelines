create_infection_compl_mat <- function (convolution_matrices,
                                        jurisdictions) {

    # with_named_dates <- lapply(convolution_matrices, function(x) {
    #     row.names(x) <- as.character(infection_days)
    #     colnames(x) <- as.character(infection_days)
    #     x
    # })
    with_named_jurisdictions <- setNames(convolution_matrices, #with_named_dates,
                                         jurisdictions)
    inf_compl_prob_mat <- do.call(cbind, lapply(with_named_jurisdictions, colSums))
    inf_compl_prob_mat
}
