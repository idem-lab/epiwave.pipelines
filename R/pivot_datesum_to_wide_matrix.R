pivot_datesum_to_wide_matrix <- function (local_summary,
                                          test_type,
                                          jurisdictions,
                                          target_dates) {

    subset <- local_summary |>
        tidyr::pivot_wider(id_cols = date_confirmation,
                           names_from = state,
                           values_from = test_type,
                           values_fn = sum,
                           values_fill = 0) |>
        data.frame()

    nonzero_dates <- subset[rowSums(subset[, -1]) != 0,]

    hh <- data.frame(
        date_confirmation = seq(min(nonzero_dates$date_confirmation),
                                max(nonzero_dates$date_confirmation), by = "days"))
    all_dates <- merge(subset, hh,
                       by = 'date_confirmation',
                       all.x = FALSE, all.y = TRUE)

    all_dates[is.na(all_dates)] <- 0

    row.names(all_dates) <- all_dates$date_confirmation
    mat <- as.matrix(all_dates[, -1])

    mat_subset <- mat[rownames(mat) %in% target_dates,]
    mat_ordered <- mat_subset[, jurisdictions]

    return(mat_ordered)
}

