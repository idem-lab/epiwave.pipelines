array_test <- array(1:120, dim = c(4, 5, 6))
dimnames(array_test) <- list(as.Date(c('2023-11-02', '2023-11-03',
                               '2023-11-04', '2023-11-05')),
                             c(paste0('draws', 1:5)),
                             NULL)

is(array_test, 'matrix')


t_test <- aperm(array_test, c(2, 1, 3))

matrix_test <- matrix(1:20, nrow = 4, ncol = 5)
dimnames(matrix_test) <- list(c('a', 'b', 'c', 'd'),
                              c(paste0('draws', 1:5)))

is(matrix_test, 'array')
