# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

# Set target options:
# packages that your targets need to run
tar_option_set(
    packages = c('greta', 'greta.gp', 'tensorflow', 'greta.dynamics'),
    format = "rds" # default storage format
)

# module <- greta::.internals$utils$misc$module

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# future::plan(future.callr::callr)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

values <- tibble::tibble(test = c('PCR', 'RAT'))

# Replace the target list below with your own:
upstream_targets <- list(
    tar_target(
        target_dates,
        as.character(
            seq.Date(as.Date("2022-05-18"),
                     as.Date("2022-10-18"),
                     by = "day"))
    ),
    tar_target(
        jurisdictions,
        colnames(PCR_matrix)
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
        nishiura_samples,
        readr::read_csv(
            file = "data/nishiura_samples.csv",
            col_types = readr::cols(param1 = readr::col_double(),
                                    param2 = readr::col_double()))
    ),
    tar_target(
        generation_interval_distribution,
        make_generation_interval_cdf(
            nishiura_samples)
    )
)
mapped <- tar_map(
    unlist = FALSE, # Return a nested list from tar_map()
    values = values,
    tar_target(
        delay_dist_mat,
        prepare_delay_input(
            target_dates, jurisdictions,
            constant_delay = list(ECDF_delay_constant_PCR))
    ),
    tar_target(
        infection_days,
        calculate_days_infection(delay_dist_mat)
    ),
    tar_target(
        notification_delay_distribution,
        extend_delay_mat(
            delay_dist_mat,
            infection_days,
            incubation_period)
    ),
    tar_target(
        timevarying_CAR,
        prepare_ascertainment_input(
            infection_days, jurisdictions,
            constant_ascertainment = 1,
            test_type = test)
    ),
    tar_target(
        notification_model_objects,
        create_model_notification_data(
            infections_timeseries = infection_model_objects$infections_timeseries,
            timevarying_delay_dist = notification_delay_distribution,
            timevarying_proportion = timevarying_CAR,
            timeseries_data = summary_matrix)
    ),
    tar_target(
        infection_completion_prob_mat,
        create_infection_compl_mat(
            notification_model_objects$convolution_matrices,
            jurisdictions)
    )
)
downstream_targets <- list(
    tar_combine(
        delay_list,
        mapped[['notification_delay_distribution']],
        command = list(!!!.x)
    ),
    tar_target(
        infection_days,
        calculate_days_infection(delay_list)
    ),
    tar_target(
        n_days_infection,
        length(infection_days)
    ),
    tar_target(
        infection_model_objects,
        create_infection_timeseries(
            n_jurisdictions,
            n_days_infection)
    ),
    tar_target(
        reff_model_objects,
        estimate_reff(
            infections_timeseries = infection_model_objects$infections_timeseries,
            generation_interval_mass_fxns = generation_interval_distribution)
    ),
    tar_combine(
        combined_model_objects,
        mapped[['notification_model_objects']],
        command = c(!!!.x, infection_model_objects, reff_model_objects)
    ),
    tar_target(
        fit,
        fit_model(combined_model_objects,
                  n_chains = 4,
                  max_convergence_tries = 1,
                  warmup = 500,
                  init_n_samples = 1000,
                  iterations_per_step = 1000)
    )
)
list(upstream_targets, mapped, downstream_targets)
