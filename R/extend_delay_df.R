
extend_delay_df <- function(delay_matrix) {

    #infection burnin in days, defined by the max notification delay at the start
    n_burnin <- max(which(delay_matrix[[1]](0:28) != 0)
    ) # how long is longest possible delay
    #infection burnin out days, defined by the max notification delay at the end
    n_burnout <- max(which(delay_matrix[[nrow(delay_matrix)]](0:28) != 0)
    ) # how long is longest possible delay

    delay_matrix <- rbind(
        do.call("rbind",
                replicate(n_burnin,
                          delay_matrix[1,],
                          simplify = FALSE)
        ),
        delay_matrix
    )

    delay_matrix <- rbind(
        delay_matrix,
        do.call("rbind",
                replicate(n_burnout,
                          delay_matrix[nrow(delay_matrix),],
                          simplify = FALSE)
        )
    )

    delay_matrix
}
