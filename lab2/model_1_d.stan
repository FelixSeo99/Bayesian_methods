
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  
//  real y_new;
}

model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  y ~ normal(alpha + beta * x, sigma);
  
  //y_new ~ normal(alpha + x_new * beta, sigma);
}

generated quantities {
  real sresid_apr[N];
  real cpo_apr[N];
  for (i in 1:N) {
    sresid_apr[i] = (y[i]- (alpha + beta * x[i])) / sigma;
    cpo_apr[i] = 1 / (sigma*2.506628) * exp(-0.5*((y[i] - (alpha + beta * x[i]))^2));
  }
}
  

