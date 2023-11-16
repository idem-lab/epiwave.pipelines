library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

R.utils::sourceDirectory('R/')

module <- greta::.internals$utils$misc$module

linelist_file <- "data-raw/processed_linelist_20231102.rds"
linelist <- readRDS(linelist_file)

local_summary <- summarise_linelist(linelist,
                                    import_status_option = 'local')

# make target dates
target_dates <- as.character(
                  seq.Date(as.Date("2023-07-01"),
                         as.Date("2023-11-01"),
                         by = "day"))
jurisdictions <- unique(local_summary$state)



# generic workflow from here
PCR_matrix <- pivot_datesum_to_wide_matrix(
    local_summary, 'PCR', target_dates, jurisdictions)
RAT_matrix <- pivot_datesum_to_wide_matrix(
    local_summary, 'RAT', target_dates, jurisdictions)

# get test type proportion matrices
PCR_prop_matrix <- pivot_test_type_prop_to_wide_matrix(
  local_summary, 'prop_PCR', target_dates)
RAT_prop_matrix <- pivot_test_type_prop_to_wide_matrix(
  local_summary, 'prop_RAT', target_dates)


PCR_valid_mat <- PCR_prop_matrix > 0
# make a valid check matrix for switching off RAT dates
RAT_valid_mat <- make_RAT_validity_matrix(RAT_matrix)
# condition date validity on at least 1 case being reported
RAT_valid_mat <- RAT_valid_mat & RAT_prop_matrix > 0





n_jurisdictions <- length(jurisdictions)

incubation_period_distribution <- make_incubation_period_cdf(strain = "Omicron")

## timevarying delays
# delay_dist_mat_PCR <- prepare_delay_input(
#         target_dates, jurisdictions,
#         delay_data = linelist[linelist$test_type == 'PCR',])
#
# delay_dist_mat_RAT <- prepare_delay_input(
#     target_dates, jurisdictions,
#     delay_data = linelist[linelist$test_type == 'RAT',])

## static delays
ECDF_delay_constant_PCR <- readRDS("data/ECDF_delay_constant_PCR.rds")
delay_dist_mat_PCR <- prepare_delay_input(
    target_dates, jurisdictions,
    constant_delay = list(ECDF_delay_constant_PCR))

ECDF_delay_constant_RAT <- readRDS("data/ECDF_delay_constant_RAT.rds")
delay_dist_mat_RAT <- prepare_delay_input(
    target_dates, jurisdictions,
    constant_delay = list(ECDF_delay_constant_RAT))

# combine the incubation and notification delays
# this function only works if there are no "null" in the delay_dist_mats.
# therefore revert_to_national must be true

## delay distribution
PCR_infection_days <- calculate_days_infection(delay_dist_mat_PCR)
PCR_notification_delay_distribution <- extend_delay_mat(
    delay_dist_mat_PCR,
    PCR_infection_days,
    incubation_period)

RAT_infection_days <- calculate_days_infection(delay_dist_mat_RAT)
RAT_notification_delay_distribution <- extend_delay_mat(
    delay_dist_mat_RAT,
    RAT_infection_days,
    incubation_period)

# CAR is non optional
timevarying_CAR_PCR <- prepare_ascertainment_input(
  PCR_infection_days, jurisdictions,
  constant_ascertainment = 1,
  test_type = "PCR")

timevarying_CAR_RAT <- prepare_ascertainment_input(
  RAT_infection_days, jurisdictions,
  constant_ascertainment = 1,
  test_type = "RAT")

# dow is optional
dow_correction_PCR <- create_dow_correction_objects(
    PCR_infection_days,
    n_jurisdictions,
    dataID = 'pcr')

dow_correction_RAT <- create_dow_correction_objects(
    RAT_infection_days,
    n_jurisdictions,
    dataID = 'rat')

# combine proportion objects
timevarying_proportion_PCR <- prepare_proportion_correction(
    target_dates,
    PCR_infection_days,
    timevarying_CAR_PCR,
    PCR_prop_matrix,
    dow_correction_PCR$pcr_dow_correction)

timevarying_proportion_RAT <- prepare_proportion_correction(
    target_dates,
    RAT_infection_days,
    timevarying_CAR_RAT,
    RAT_prop_matrix,
    dow_correction_RAT$rat_dow_correction)

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
    observable_infection_dates = PCR_infection_days,
    timevarying_delay_dist = PCR_notification_delay_distribution,
    timevarying_proportion = timevarying_proportion_PCR$timevarying_proportion,
    observed_data = PCR_matrix,
    valid_mat = PCR_valid_mat,
    dataID = 'pcr')

RAT_notification_model_objects <- create_model_notification_data(
    infections_timeseries = infection_model_objects$infections_timeseries,
    full_infection_dates = days_infection,
    observable_infection_dates = RAT_infection_days,
    timevarying_delay_dist = RAT_notification_delay_distribution,
    timevarying_proportion = timevarying_proportion_RAT$timevarying_proportion,
    observed_data = RAT_matrix,
    valid_mat = RAT_valid_mat,
    dataID = 'rat')

# priors for the parameters of the lognormal distribution over the serial
#interval from Nishiura et al., as stored in the EpiNow source code
nishiura_samples <- readr::read_csv(
    file = "data/nishiura_samples.csv",
    col_types = readr::cols(param1 = readr::col_double(),
                            param2 = readr::col_double()))

generation_interval_distribution <- make_generation_interval_density(
    nishiura_samples)

reff_model_objects <- estimate_reff(
    infections_timeseries = infection_model_objects$infections_timeseries,
    generation_interval_mass_fxns = generation_interval_distribution)

combined_model_objects <- c(timevarying_proportion_PCR,
                            timevarying_proportion_RAT,
                            infection_model_objects,
                            PCR_notification_model_objects,
                            RAT_notification_model_objects,
                            reff_model_objects)

m <- model(combined_model_objects$infections_timeseries,
           combined_model_objects$reff)
plot(m)

fit <- fit_model(model = m,
                 n_chains = 4,
                 max_convergence_tries = 1,
                 warmup = 1000,
                 init_n_samples = 1000,
                 iterations_per_step = 1000)



###=== IN DEV

# infection completion probability matrices
PCR_infection_completion_prob_mat <- create_infection_compl_mat(
  observable_infection_dates = PCR_infection_days,
  target_dates = target_dates,
  jurisdictions = jurisdictions,
  timevarying_delay_dist_ext = PCR_notification_delay_distribution)

RAT_infection_completion_prob_mat <- create_infection_compl_mat(
  observable_infection_dates = RAT_infection_days,
  target_dates = target_dates,
  jurisdictions = jurisdictions,
  timevarying_delay_dist_ext = RAT_notification_delay_distribution)

# #check convergence
# coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]


case_sims_RAT <- calculate(combined_model_objects$rat_observed_data_array,
                           values = fit,
                           nsim = 1000)

plot_timeseries_sims(case_sims_RAT[[1]],
                     type = "notification",
                     dates = as.Date(rownames(RAT_matrix)),
                     states = colnames(RAT_matrix),
                     valid_mat = RAT_valid_mat,
                     start_date = as.Date("2023-05-01"),
                     dim = "1",
                     case_validation_data = local_summary |>
                       dplyr::rename("date" = date_confirmation,
                                     "count" = RAT))

case_sims_PCR <- calculate(combined_model_objects$pcr_observed_data_array,
                           values = fit,
                           nsim = 1000)

plot_timeseries_sims(case_sims_PCR[[1]],
                     type = "notification",
                     dates = as.Date(rownames(PCR_matrix)),
                     states = colnames(PCR_matrix),
                     valid_mat = PCR_valid_mat,
                     start_date = as.Date(rownames(PCR_matrix)[1]),
                     dim = "1",
                     case_validation_data = local_summary |>
                       dplyr::rename("date" = date_confirmation,
                                     "count" = PCR))

infection_sims <- calculate(combined_model_objects$infections_timeseries,
                            values = fit,
                            nsim = 1000)

plot_timeseries_sims(infection_sims[[1]],
                     type = "infection",
                     dates = days_infection,
                     start_date = as.Date(rownames(PCR_matrix)[1]),
                     states = jurisdictions,
                     dim_sim = "2",
                     case_validation_data = local_summary |>
                       dplyr::rename("date" = date_confirmation,
                                     "count" = total))

reff_sims <- calculate(combined_model_objects$reff,
                       values = fit,
                       nsim = 1000)

plot_timeseries_sims(reff_sims[[1]],
                     type = "reff",
                     dates = days_infection,
                     start_date = as.Date(rownames(PCR_matrix)[1]),
                     end_date = as.Date(rownames(PCR_matrix)[nrow(PCR_matrix)]),
                     states = jurisdictions, dim_sim = "2")


calculate(combined_model_objects$gp_lengthscale,
          values = fit,
          nsim = 10)
calculate(combined_model_objects$gp_variance,
          values = fit,
          nsim = 10)
