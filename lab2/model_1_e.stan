data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real alpha;
  real beta_1;
  real beta_2;
  real<lower=0> sigma;

}

model {
  alpha ~ normal(0, 10);
  beta_1 ~ normal(0, 10);
  beta_2 ~ normal(0, 10);
  sigma ~ normal(0, 10);
  y ~ normal(alpha + beta_1 * x + beta_2 * x^2, sigma);
  
}

  