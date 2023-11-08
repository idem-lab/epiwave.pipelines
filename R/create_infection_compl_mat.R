create_infection_compl_mat <- function (observable_infection_dates,
                                        target_dates,
                                        jurisdictions,
                                        timevarying_delay_dist_ext) {

    # get length of infection timeseries
    n_days_infection_observable <- length(observable_infection_dates)

    # make index for target days within the infection days
    # basically this backs out extra left and right days
    obs_data_idx <- which(observable_infection_dates %in% as.Date(target_dates))

    # build convolution matrix
    convolution_matrices <- lapply(1:n_jurisdictions, function(x)
        get_convolution_matrix(timevarying_delay_dist_ext[, x],
                               n_days_infection_observable))

    # subset convolution matrix to target days
    convolution_matrices <- lapply(convolution_matrices,
                                   function(x){
                                       x[obs_data_idx,obs_data_idx]
                                   })
    # set jurisdiction names
    with_named_jurisdictions <- setNames(convolution_matrices, #with_named_dates,
                                         jurisdictions)
    # summarise completion probability by state
    inf_compl_prob_mat <- do.call(cbind, lapply(with_named_jurisdictions, colSums))

    # put in date names
    rownames(inf_compl_prob_mat) <- target_dates
    # retun date by state object
    inf_compl_prob_mat
}
