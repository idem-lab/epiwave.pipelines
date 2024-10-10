n_states <- 8
n_days <- 100

# fake parameters
infections_eps <- matrix(rnorm(n_days * n_states, 0, 0.1), n_days)
infections_mn <- 5 + 1 * sin(seq_len(n_days) / 14)
infections <- exp(sweep(infections_eps, 1, infections_mn, FUN = "+"))



not_to_infection_delays <- XXXX
convolution_matrix_cases <- get_convolution_matrix(
    not_to_infection_delays, nrow(infections))
ave_cases <- convolution_matrix_cases %*% infections[,1] * CAR
cases <- rpois(n_days * n_states, ave_cases) # reshape


nishiura_samples <- readr::read_csv(
    file = "data/nishiura_samples.csv",
    col_types = readr::cols(param1 = readr::col_double(),
                            param2 = readr::col_double()))
generation_interval_distribution <- make_generation_interval_density(
    nishiura_samples)
convolution_matrix <- get_convolution_matrix(
    generation_interval_distribution, nrow(infections))
infectiousness <- convolution_matrix %*% infections[,1]

reff <- infections[,1] / infectiousness

plot(1:100, infections[,1])
plot(4:100, reff[4:100])
plot(1:100, infectiousness[,1])

# then can use cases and delays in new code to est. reff
