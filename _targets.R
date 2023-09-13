# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint
# source('scripts/data_preprocessing.R')
# Set target options:
# packages that your targets need to run
tar_option_set(
    packages = c('greta', 'greta.gp', 'tensorflow', 'greta.dynamics'),
    format = "rds" # default storage format
    # Set other options as needed.
)

module <- greta::.internals$utils$misc$module

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(
    tar_target(
        notification_matrix,
        readRDS('tests/sim_case_data.RData')
    ),
    tar_target(
        hospitalisation_matrix,
        # notification_matrix + 10 * .5
        readRDS('tests/sim_hospital_data.rds')
    ), # connect this to the others
    tar_target(
        jurisdictions,
        colnames(notification_matrix)
    ),
    tar_target(
        n_jurisdictions,
        length(jurisdictions)
    ),
    tar_target(
        incubation_period_distribution,
        make_incubation_period_cdf(strain = "Omicron")
    ),
    tar_target(
        notification_delay_distribution,
        create_timevarying_delay_distribution(
            notification_matrix,
            incubation_period_distribution)
    ),
    # creating delay data: see "rolling_delays.R" in old reff code
    tar_target(
        hospitalisation_delay_distribution,
        readRDS('tests/hospitalisation_delay_distribution.rds')
    ),
    tar_target(
        n_days_infection,
        calculate_n_days_infection(
            list(notification_delay_distribution,
                 hospitalisation_delay_distribution))
    ),
    tar_target(
        timevarying_CAR, # link with CAR matrix in BDSS repo
        prepare_ascertainment_input(
            assume_constant_ascertainment = TRUE,
            constant_ascertainment_fraction = 1)
    ),
    tar_target(
        timevarying_proportion_hospitalised,
        prepare_ascertainment_input(
            assume_constant_ascertainment = TRUE,
            constant_ascertainment_fraction = 1)
    ),
    tar_target(
        infection_model_objects,
        create_infection_timeseries(
            n_jurisdictions,
            n_days_infection)
    ),
    tar_target(
        notification_model_objects,
        create_negbin_model_notification_data(
            infection_model_objects,
            timeseries_data = notification_matrix,
            timevarying_delay_distribution = notification_delay_distribution,
            timevarying_proporion = timevarying_CAR,
            n_days_infection,
            n_jurisdictions)
    ),
    tar_target(
        hospitalisation_model_objects,
        create_poisson_model_hospitalisation_data(
            infection_model_objects,
            timeseries_data = hospitalisation_matrix,
            timevarying_delay_distribution = hospitalisation_delay_distribution,
            timevarying_proporion = timevarying_proportion_hospitalised,
            n_days_infection,
            n_jurisdictions)
    ),
    tar_target(
        combined_model_objects,
        c(infection_model_objects,
          notification_model_objects,
          hospitalisation_model_objects)
    ),
    tar_target(
        fit,
        fit_model(combined_model_objects,
                  n_chains = 4,
                  max_convergence_tries = 2,
                  warmup = 500,
                  init_n_samples = 1000,
                  iterations_per_step = 1000) # this doesn't feel like it needs to be user defined?
    )
)
