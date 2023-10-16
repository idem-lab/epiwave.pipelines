library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

R.utils::sourceDirectory('R/')

module <- greta::.internals$utils$misc$module

linelist_file <- "data-raw/linelist_processed_2023_09_21.rds"
linelist <- readRDS(linelist_file)

local_summary <- summarise_linelist(linelist,
                                    import_status_option = 'local')

#make target dates for end of RAT dates
target_dates <- as.character(seq.Date(as.Date("2022-03-01"),as.Date("2022-08-01"),by = "day"))


PCR_matrix <- pivot_datesum_to_wide_matrix(local_summary, 'PCR')
PCR_matrix <- PCR_matrix[rownames(PCR_matrix) %in% target_dates,]
RAT_matrix <- pivot_datesum_to_wide_matrix(local_summary, 'RAT')
RAT_matrix <- RAT_matrix[rownames(RAT_matrix) %in% target_dates,]
## ensure all have all jurisdictions even if some test types only have
jurisdictions <- unique(linelist$state)

# set state names
# if (is.null(jurisdiction_names)) {
#     jurisdictions <- colnames(notification_matrix)
# } else if (length(jurisdiction_names) != ncol(notification_matrix)) {
#     stop("Error: supplied jurisdiction names has different length from number of columns in notification matrix")
# } else {
#     jurisdictions <- jurisdiction_names
# }

n_jurisdictions <- length(jurisdictions)

incubation_period_distribution <- make_incubation_period_cdf(strain = "Omicron")

# testing new delay function
# make this separate for pcr and rat
# delay_dist_mat_PCR <- estimate_delays(
#     linelist[linelist$test_type == 'PCR',],
#     jurisdictions)
# targets::tar_load(delay_dist_mat_PCR)
# delay_dist_mat_RAT <- estimate_delays(
#     linelist[linelist$test_type == 'RAT',],
#     jurisdictions)
#targets::tar_load(delay_dist_mat_RAT)

ECDF_delay_constant_PCR <- readRDS("data/ECDF_delay_constant_PCR.rds")
delay_dist_mat_PCR <- matrix(list(ECDF_delay_constant_PCR),
                             nrow = nrow(PCR_matrix),
                             ncol = n_jurisdictions,
                             dimnames = list(rownames(PCR_matrix),
                                             jurisdictions))
ECDF_delay_constant_RAT <- readRDS("data/ECDF_delay_constant_RAT.rds")
delay_dist_mat_RAT <- matrix(list(ECDF_delay_constant_RAT),
                             nrow = nrow(RAT_matrix),
                             ncol = n_jurisdictions,
                             dimnames = list(rownames(RAT_matrix),
                                             jurisdictions))

# apply construct_delays to each cell
# combine the incubation and notification delays
# this function only works if there are no "null" in the delay_dist_mats.
# therefore revert_to_national must be true

# delay distribution
PCR_notification_delay_distribution <- apply(
    delay_dist_mat_PCR,
    c(1,2), construct_delays,
    ecdf2 = incubation_period_distribution,
    output = "probability",
    stefun_output = TRUE)

PCR_infection_days <- calculate_days_infection(PCR_notification_delay_distribution)

PCR_notification_delay_distribution_ext <- extend_delay_df(PCR_notification_delay_distribution)
rownames(PCR_notification_delay_distribution_ext) <- as.character(PCR_infection_days)

RAT_notification_delay_distribution <- apply(
    delay_dist_mat_RAT,
    c(1,2), construct_delays,
    ecdf2 = incubation_period_distribution,
    output = "probability",
    stefun_output = TRUE)

RAT_infection_days <- calculate_days_infection(RAT_notification_delay_distribution)

RAT_notification_delay_distribution_ext <- extend_delay_df(RAT_notification_delay_distribution)
rownames(RAT_notification_delay_distribution_ext) <- as.character(RAT_infection_days)


timevarying_CAR_PCR <- prepare_ascertainment_input(
  PCR_infection_days, jurisdictions,constant_ascertainment = 1,
    test_type = "PCR")

timevarying_CAR_RAT <- prepare_ascertainment_input(
  RAT_infection_days, jurisdictions,
  constant_ascertainment =
        1,
    test_type = "RAT")

############# above objects are all created based on input data, and does not
############# require the creation of an infection timeseries

days_infection <- seq.Date(min(PCR_infection_days,
                               RAT_infection_days),
                           max(PCR_infection_days,
                               RAT_infection_days),
                           by = "day")

n_days_infection <- length(days_infection)

infection_model_objects <- create_infection_timeseries(
    n_jurisdictions,
    n_days_infection)

PCR_notification_model_objects <- create_model_notification_data(
    infections_timeseries = infection_model_objects$infections_timeseries,
    full_infection_dates = days_infection,
    observed_infection_dates = PCR_infection_days,
    timevarying_delay_dist = PCR_notification_delay_distribution_ext,
    timevarying_proportion = timevarying_CAR_PCR,
    timeseries_data = PCR_matrix)

RAT_notification_model_objects <- create_model_notification_data(
    infections_timeseries = infection_model_objects$infections_timeseries,
    full_infection_dates = days_infection,
    observed_infection_dates = RAT_infection_days,
    timevarying_delay_dist = RAT_notification_delay_distribution_ext,
    timevarying_proportion = timevarying_CAR_RAT,
    timeseries_data = RAT_matrix)

# priors for the parameters of the lognormal distribution over the serial
#interval from Nishiura et al., as stored in the EpiNow source code
nishiura_samples <- readr::read_csv(
    file = "data/nishiura_samples.csv",
    col_types = readr::cols(param1 = readr::col_double(),
                            param2 = readr::col_double()))

generation_interval_distribution <- make_generation_interval_cdf(
    nishiura_samples)

generation_interval_distribution <- construct_delays(generation_interval_distribution,
                                                     output = "probability",
                                                     stefun_output = TRUE)

reff_model_objects <- estimate_reff(
    infections_timeseries = infection_model_objects$infections_timeseries,
    generation_interval_mass_fxns = generation_interval_distribution)

combined_model_objects <- c(infection_model_objects,
                            PCR_notification_model_objects,
                            RAT_notification_model_objects,
                            reff_model_objects)

fit <- fit_model(combined_model_objects,
                 n_chains = 2,
                 max_convergence_tries = 1,
                 warmup = 400,
                 init_n_samples = 1000,
                 iterations_per_step = 1000) # this doesn't feel like it needs to be user defined?

# reff <- calculate(fit)

###=== IN DEV

# infection completion probability matrices
PCR_infection_completion_prob_mat <- create_infection_compl_mat(
    PCR_notification_model_objects$convolution_matrices,
    jurisdictions)

RAT_infection_completion_prob_mat <- create_infection_compl_mat(
    RAT_notification_model_objects$convolution_matrices,
    jurisdictions)

# #check convergence
# coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]


case_sims <- calculate(combined_model_objects$timeseries_data_array,
                       values = fit,
                       nsim = 1000)

plot_timeseries_sims(case_sims[[1]],
                     type = "notification",
                     dates = seq.Date(as.Date("2022-03-01"),as.Date("2022-08-01"),by = "day"),
                     states = colnames(PCR_matrix))



calculate(combined_model_objects$gp_lengthscale,
          values = fit,
          nsim = 10)
calculate(combined_model_objects$gp_variance,
          values = fit,
          nsim = 10)
