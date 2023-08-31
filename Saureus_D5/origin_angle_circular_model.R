# Rscript to perform circular regression on origin-angle data

library(rstan)
library(tidyverse)

# Set options for Rstan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set working directory (modify as needed)
setwd('### SET YOUR WORKING DIRECTORY HERE ###')

# Import data and convert to the correct format
data_table <- read_tsv('tsb_origin_angle_scaled.txt')

dat <- list(
  'N' = data_table[['N']][1],
  'D' = data_table[['D']],
  'A' = data_table[['A']]
) 

# Specify and save the stan model
model <- 'origin_angle.stan'
write("// Simple linear model relating origin distance to angle
// See: https://stats.stackexchange.com/questions/95085/multiple-regression-in-directional-circular-statistics
// Also for reference without intercept: Jammalamadaka, S. R. and SenGupta, A. (2001). Topics in Circular Statistics. World Scientific, Singapore (Chapter 8).
data {
  int <lower=0> N; // Number of genes
  vector <lower=-pi(), upper=pi()> [N] D; // Origin distance
  vector <lower=-pi(), upper=pi()> [N] A; // Gene angle
}

parameters {
  // See manuscript for full equation
  // Gradient
  real <lower=0> k;
  // Intercept
  real b1;
  real b2;
  // Noise
  real <lower=0> kappa;
}

model {
  vector [N] sin_A;
  vector [N] cos_A;
  vector [N] mu;
  k ~ lognormal(0, 0.5);
  b1 ~ normal(0, 0.5);
  b2 ~ normal(0, 0.5);
  kappa ~ lognormal(0, 1.0);
  
  cos_A = b1 * cos(k*D) - b2 * sin(k*D);
  sin_A = b2 * cos(k*D) + b1 * sin(k*D);
  for (n in 1:N) {
    mu[n] = atan2(sin_A[n], cos_A[n]);
  }
  
  A ~ von_mises(mu, kappa);
}

", model)

# Compile the stan model, checking it works
stanc(model)

# Fit the stan model
fit <- stan(model, data=dat, chains=8, iter = 2000, seed=123)

# Note - chains often give different fits, presumably due to 
# local minima. It is easy to tell which are the correct fits
# based on the log posterior estimate "lp__". Pick the chains
# that fit to the highest "lp__" values.

# Plot for each chain:
traceplot(fit, pars='lp__')

# Plot which chains cross a threshold (you can set this manually)
hist(c(1:8000)[rstan::extract(fit, 'lp__')[[1]] > 500])

# Specify the intervals matching to chains with a good fit
intervals <- c(1:2000, 4001:8000)

# Extract the parameters (note if you have filtered the chains correctly,
# these should be unimodal in a histogram).
k <- rstan::extract(fit, 'k')[[1]][intervals]
b1 <- rstan::extract(fit, 'b1')[[1]][intervals]
b2 <- rstan::extract(fit, 'b2')[[1]][intervals]

# Convert these into values for prediction
cos_A <- apply(cos(k %*% t(dat$D)), 2, function(x) x * b1) +
  apply(sin(k %*% t(dat$D)), 2, function(x) (x * -b2))
sin_A <- apply(cos(k %*% t(dat$D)), 2, function(x) x * b2) +
  apply(sin(k %*% t(dat$D)), 2, function(x) (x * b1))

A_pred <- atan2(sin_A, cos_A)

# Determine mean prediction (this calculates the angular mean)
angle_mean <- function(theta) {
  return(atan2(mean(sin(theta)), mean(cos(theta))))
}


angle_pred <- apply(A_pred, 2, angle_mean)

# Get the mean k value
mean(k) # 1.400709

# Now save the outputs
write.table(data.frame(gene=data_table$...1, angle_pred=angle_pred), 'tsb_predicted_angle.txt', sep="\t", quote=F, row.names=F)

# Save model
saveRDS(fit, 'tsb_stanfit.rds')

# Predict the gene angle of a gene at position = 0
origin_position <- -pi # Because I've scaled from -pi to pi

cos_A <- apply(cos(k %*% t(origin_position)), 2, function(x) x * b1) +
  apply(sin(k %*% t(origin_position)), 2, function(x) (x * -b2))
sin_A <- apply(cos(k %*% t(origin_position)), 2, function(x) x * b2) +
  apply(sin(k %*% t(origin_position)), 2, function(x) (x * b1))

A_pred <- atan2(sin_A, cos_A)

angle_mean(A_pred) # Gives 2.908317

# Next, predict when the terminus is reached:
origin_position <- pi # Because I've scaled from -pi to pi

cos_A <- apply(cos(k %*% t(origin_position)), 2, function(x) x * b1) +
  apply(sin(k %*% t(origin_position)), 2, function(x) (x * -b2))
sin_A <- apply(cos(k %*% t(origin_position)), 2, function(x) x * b2) +
  apply(sin(k %*% t(origin_position)), 2, function(x) (x * b1))

A_pred <- atan2(sin_A, cos_A)

angle_mean(A_pred) # Gives -0.8571364

# Polymerase speed
chrom_length <- 2872769
td = 24.9 # Doubling time in min

((chrom_length / 2) / mean(k)) / (td * 60) # 686.3918 bp/s
