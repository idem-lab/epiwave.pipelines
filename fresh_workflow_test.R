library(tidyverse)
library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

R.utils::sourceDirectory('R/')

module <- greta::.internals$utils$misc$module

#test with sim data
notification_matrix <- readRDS("tests/sim_case_data.RData")
# add ten days, multiple by .5
jurisdictions <- colnames(notification_matrix)

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

# shoudl be same dims as not. matrix but not. matrix not necessarily as input
notification_delay_distribution <- create_timevarying_delay_distribution(
    notification_matrix,
    incubation_period_distribution)

n_days_infection <- calculate_n_days_infection(
    list(notification_delay_distribution))

timevarying_CAR <- prepare_ascertainment_input(
    assume_constant_ascertainment = TRUE,
    constant_ascertainment_fraction = 1)

infection_model_objects <- create_infection_timeseries(
    n_jurisdictions,
    n_days_infection)

# debug
notification_model_objects <- create_negbin_model_notification_data(
    gp_model_objects,
    timeseries_data = notification_matrix,
    timevarying_delay_distribution = notification_delay_distribution,
    timevarying_proporion = timevarying_CAR,
    n_days_infection,
    n_jurisdictions)

combined_model_objects <- c(infection_model_objects, notification_model_objects)

fit <- fit_model(combined_model_objects,
                 n_chains = 4,
                 max_convergence_tries = 1,
                 warmup = 500,
                 init_n_samples = 1000,
                 iterations_per_step = 1000) # this doesn't feel like it needs to be user defined?
