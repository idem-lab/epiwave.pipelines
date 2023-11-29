# header
library(greta.gp)
library(greta)
library(tensorflow)
library(greta.dynamics)
R.utils::sourceDirectory('R/')
module <- greta::.internals$utils$misc$module

source('sitassess_processing.R')

# user input
pcr_dat <- clean_pcr
constant_delay_data <- readRDS("data/ECDF_delay_constant_PCR.rds")
timevarying_delay_data <- NULL
constant_ascertainment <- 1
data_id <- 'pcr'

dates <- unique(pcr_dat$date)
jurisdictions <- unique(pcr_dat$jurisdiction)
ECDF_delay_constant_PCR <- readRDS("data/ECDF_delay_constant_PCR.rds")
delay_dist_mat_PCR <- prepare_delay_input(
    dates, jurisdictions,
    constant_delay = list(ECDF_delay_constant_PCR))
infection_days <- calculate_days_infection(delay_dist_mat_PCR)

delay_input <- matrix(list(constant_delay_data),
                      nrow = 212,
                      ncol = 8,
                      dimnames = list(as.character(infection_days),
                                      jurisdictions))
# infection_days <- seq.Date(min(as.Date(target_dates))-40,
#                            max(as.Date(target_dates))+40,
#                            by = 'day')

n_jurisdictions <- length(unique(pcr_dat$jurisdiction))

long_unique <- expand.grid(date = infection_days,
                           jurisdiction = unique(pcr_dat$jurisdiction))

pcr_car <- cbind(
    long_unique,
    tibble::tibble(
        prop = rep(constant_ascertainment, nrow = nrow(clean_pcr))
    )
)

notification_delay_dist <- cbind(
    long_unique,
    tibble::tibble(
        delay_fxn = rep(list(constant_delay_data),
                        nrow = length(infection_days))
    )
)

# combine the incubation and notification delays
# this function only works if there are no "null" in the delay_dist_mats.
# therefore revert_to_national must be true
incubation_period_distribution <- lowerGPreff::make_incubation_period_cdf(
    strain = "Omicron")

## delay distribution
infection_delay_distribution <- lowerGPreff::extend_delay_data(
    delay_input,
    infection_days,
    incubation_period_distribution)

# dow is optional
dow_model <- lowerGPreff::create_dow_correction_arrays(
    n_jurisdictions,
    data_id = data_id)


# timevarying_CAR <- get_CAR_from_surveys(
#     dates, jurisdictions, ascertainment_estimate,
#     test_type = 'PCR')

# # combine proportion objects
# timevarying_proportion <- prepare_proportion_correction(
#     dates,
#     input_specific_infection_days,
#     timevarying_ascertainment = timevarying_CAR,
#     # case_type_proportion = case_data_diagnostic_output$PCR_prop_matrix,
#     dow_correction = dow_correction$pcr_dow_correction)


n_days_infection <- length(infection_days)

infection_model_objects <- lowerGPreff::create_infection_timeseries(
    n_days_infection,
    n_jurisdictions,
    effect_type = 'growth_rate')

PCR_observation_model_objects <- create_observation_model(
    # put through the whole infection_model_objects instead?
    infection_timeseries = infection_model_objects$infection_timeseries,
    infection_delay_distribution = infection_delay_distribution,
    proportion_observed = pcr_car,
    count_data = pcr_dat,
    dow_model = dow_model,
    data_id = data_id)

# priors for the parameters of the lognormal distribution over the serial
#interval from Nishiura et al., as stored in the EpiNow source code
nishiura_samples <- readr::read_csv(
    file = "data/nishiura_samples.csv",
    col_types = readr::cols(param1 = readr::col_double(),
                            param2 = readr::col_double()))

generation_interval_distribution <- lowerGPreff::make_generation_interval_density(
    nishiura_samples)

reff_model_objects <- estimate_reff(
    infection_timeseries = infection_model_objects$infection_timeseries,
    generation_interval_mass_fxns = generation_interval_distribution)

combined_model_objects <- c(timevarying_proportion,
                            # timevarying_proportion_RAT,
                            infection_model_objects,
                            PCR_observation_model_objects,
                            # RAT_notification_model_objects,
                            reff_model_objects)

m <- model(combined_model_objects$infection_timeseries,
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

# get nowcast start date
nowcast_start <- as.Date(rownames(PCR_infection_completion_prob_mat)[min(which(PCR_infection_completion_prob_mat < 0.95))])

# #check convergence
# coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]


case_sims_RAT <- calculate(combined_model_objects$rat_observed_data_array,
                           values = fit,
                           nsim = 1000)

plot_timeseries_sims(case_sims_RAT[[1]],
                     type = "notification",
                     dates = as.Date(rownames(RAT_matrix)),
                     states = colnames(RAT_matrix),
                     valid_mat = case_data_diagnostic_output$RAT_valid_mat,
                     start_date = as.Date("2023-05-01"),
                     dim = "1",
                     case_validation_data = local_summary |>
                         dplyr::rename("date" = date_confirmation,
                                       "count" = RAT),
                     nowcast_start = nowcast_start)

case_sims_PCR <- calculate(combined_model_objects$pcr_observed_data_array,
                           values = fit,
                           nsim = 1000)

plot_timeseries_sims(case_sims_PCR[[1]],
                     type = "notification",
                     dates = as.Date(rownames(PCR_matrix)),
                     states = colnames(PCR_matrix),
                     valid_mat = case_data_diagnostic_output$PCR_valid_mat,
                     start_date = as.Date(rownames(PCR_matrix)[1]),
                     dim = "1",
                     case_validation_data = local_summary |>
                         dplyr::rename("date" = date_confirmation,
                                       "count" = PCR),
                     nowcast_start = nowcast_start)

infection_sims <- calculate(combined_model_objects$infection_timeseries,
                            values = fit,
                            nsim = 1000)

plot_timeseries_sims(infection_sims[[1]],
                     type = "infection",
                     dates = infection_timeseries_dates,
                     start_date = as.Date(rownames(linelist_pcr)[1]),
                     end_date = as.Date(rownames(linelist_pcr)[nrow(linelist_pcr)]),
                     states = jurisdictions,
                     dim_sim = "2",
                     infection_nowcast = FALSE)#,
                     # nowcast_start = nowcast_start)

reff_sims <- calculate(combined_model_objects$reff,
                       values = fit,
                       nsim = 1000)

plot_timeseries_sims(reff_sims[[1]],
                     type = "reff",
                     dates = infection_timeseries_dates,
                     start_date = as.Date(rownames(pcr_mat[,,'notification_data'])[1]),
                     end_date = as.Date(rownames(pcr_mat[,,'notification_data'])
                                        [nrow(pcr_mat[,,'notification_data'])]),
                     states = jurisdictions, dim_sim = "2",
                     infection_nowcast = FALSE)#,
                     # nowcast_start = nowcast_start)


calculate(combined_model_objects$gp_lengthscale,
          values = fit,
          nsim = 10)
calculate(combined_model_objects$gp_variance,
          values = fit,
          nsim = 10)
