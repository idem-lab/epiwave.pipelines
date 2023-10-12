library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

source('R/TF2.R')

# variables & priors
int <- normal(0, 10)
coef <- normal(0, 10)
sd <- cauchy(0, 3, truncation = c(0, Inf))

# linear predictor
mu <- int + coef * attitude$complaints

# observation model
distribution(attitude$rating) <- normal(mu, sd)

m <- model(int, coef, sd)

n_chains <- 4
max_convergence_tries <- 2
warmup <- 500
init_n_samples <- 1000
iterations_per_step <- 1000

# get stable inits
init <- generate_valid_inits(model = m,
                             chains = n_chains,
                             max_tries = 1000)

