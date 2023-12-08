## header
library(greta.gp)
library(greta)
library(tensorflow)
library(greta.dynamics)
R.utils::sourceDirectory('R/')

source('sitassess_processing.R')

## user input
pcr_dat <- clean_pcr
constant_delay_data <- readRDS("data/ECDF_delay_constant_PCR.rds")
timevarying_delay_data <- NULL
constant_ascertainment <- 1
data_id <- 'pcr'

infection_days <- short_term_get_infection_days(
    pcr_dat,
    constant_delay_data)
jurisdictions <- unique(pcr_dat$jurisdiction)
n_jurisdictions <- length(jurisdictions)

car <- create_long_data_for_constant_val(
    dates = infection_days,
    jurisdictions = jurisdictions,
    constant_val = 0.75,
    col_name = 'prop'
)

notification_delay_dist <- create_long_data_for_constant_val(
    dates = infection_days,
    jurisdictions = jurisdictions,
    constant_val = list(constant_delay_data),
    col_name = 'delay_fxn'
)

# priors for the parameters of the lognormal distribution over the serial
# interval from Nishiura et al., as stored in the EpiNow source code
gi_distribution_data <- readr::read_csv(
    file = "data/nishiura_samples.csv",
    col_types = readr::cols(param1 = readr::col_double(),
                            param2 = readr::col_double()))


## generic flow
# combine the incubation and notification delays
# this function only works if there are no "null" in the delay_dist_mats.
incubation_period_distribution <- lowerGPreff::make_incubation_period_cdf(
    strain = "Omicron")

# delay distribution
infection_delay_distribution <- lowerGPreff::extend_delay_data(
    notification_delay_dist,
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

PCR_observation_model_objects <- lowerGPreff::create_observation_model(
    # put through the whole infection_model_objects instead?
    infection_timeseries = infection_model_objects$infection_timeseries,
    delay_distribution = infection_delay_distribution,
    proportion_observed = car,
    count_data = pcr_dat,
    dow_model = dow_model,
    data_id = data_id)

generation_interval_distribution <- lowerGPreff::make_generation_interval_density(
    gi_distribution_data)

reff_model_objects <- lowerGPreff::estimate_reff(
    infection_timeseries = infection_model_objects$infection_timeseries,
    generation_interval_mass_fxns = generation_interval_distribution)

m <- model(infection_model_objects$infection_timeseries,
           reff_model_objects$reff)
plot(m)

fit <- fit_model(model = m,
                 n_chains = 4,
                 max_convergence_tries = 1,
                 warmup = 1000,
                 init_n_samples = 1000,
                 iterations_per_step = 1000)

#### check and revise
infection_traj <- build_infection_trajectories(
    param = infection_model_objects$infection_timeseries,
    infection_days,
    fit,
    infection_delay_distribution,
    jurisdictions)

reff_intervals <- build_reff_estimates(
    param = reff_model_objects$reff,
    infection_days,
    fit,
    jurisdictions)


#### to here ready



# get nowcast start date
nowcast_start <- Reduce(min, lapply(infections_out, function (x)
    as.Date(rownames(x)[min(which(x$inf_compl < 0.95))])
))

reff_estimates





# #check convergence
# coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]



MCMCvis::MCMCsummary(infection_sims)


plot_timeseries_sims(infection_sims[[1]],
                     type = "infection",
                     dates = infection_timeseries_dates,
                     start_date = as.Date(rownames(linelist_pcr)[1]),
                     end_date = as.Date(rownames(linelist_pcr)[nrow(linelist_pcr)]),
                     states = jurisdictions,
                     dim_sim = "2",
                     infection_nowcast = FALSE)#,
# nowcast_start = nowcast_start)

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
