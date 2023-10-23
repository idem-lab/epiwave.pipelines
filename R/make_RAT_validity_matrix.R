#' takes a RAT notification matrix, returns a logical matrix designating if RAT
#' data is available and reliable for that date and state, based on when RAT
#' reporting was officially ceased for a state. The input matrix must have state
#' as colnames and date as rownames, the date for off-switch is hard-coded for
#' each state
#'
#' @param RAT_matrix
#'
#' @return
#' @export
#'
#' @examples
make_RAT_validity_matrix <- function(RAT_matrix) {

    RAT_valid_mat <- RAT_matrix
    RAT_valid_mat[] <- TRUE
    RAT_valid_mat[as.Date(rownames(RAT_valid_mat)) > as.Date("2023-06-30"),
                  colnames(RAT_valid_mat) == "VIC"] <- FALSE
    RAT_valid_mat[as.Date(rownames(RAT_valid_mat)) > as.Date("2023-08-23"),
                  colnames(RAT_valid_mat) == "SA"] <- FALSE
    RAT_valid_mat[as.Date(rownames(RAT_valid_mat)) > as.Date("2023-09-01"),
                  colnames(RAT_valid_mat) == "QLD"] <- FALSE
    #coerce to boolean
    RAT_valid_mat <- ifelse(RAT_valid_mat == 1,TRUE, FALSE)
}
