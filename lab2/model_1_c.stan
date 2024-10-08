
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  
  real x_new;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  
//  real y_new;
}

//transformed parameters {
//  real<lower=0,upper=1> p;
  
//  p = 1 - normal_cdf(0, alpha + beta * x_new, sigma);
//}

model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  y ~ normal(alpha + beta * x, sigma);
  
  //y_new ~ normal(alpha + x_new * beta, sigma);
}

generated quantities {
  real y_new = normal_rng(alpha + beta * x_new, sigma);
}
