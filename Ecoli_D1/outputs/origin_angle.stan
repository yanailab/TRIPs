// Simple linear model relating origin distance to angle
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


