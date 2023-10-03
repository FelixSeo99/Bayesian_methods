
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  
  int<lower=0> N_new;
  real x_new;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  
  real y_new;
}

transformed parameters {
  real<lower=0,upper=1> p;
  
  p = 1 - normal_cdf(0, alpha + beta * x, sigma);
}

model {
  alpha ~ uniform(-100, 100);
  beta ~ uniform(-100, 100);
  sigma ~ uniform(-100, 100);
  y ~ normal(alpha + beta * x, sigma);
  
  y_new ~ normal(alpha + x_new * beta, sigma);
}

//generated quantities {
//  real y_rep[N] = normal_rng(alpha + beta * x, sigma);
//}
