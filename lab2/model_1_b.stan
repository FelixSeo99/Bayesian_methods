
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ uniform(-100, 100);
  beta ~ uniform(-100, 100);
  sigma ~ uniform(-100, 100);
  y ~ normal(alpha + beta * x, sigma);
}

