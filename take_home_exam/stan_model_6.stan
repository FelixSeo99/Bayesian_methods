
data {
  int<lower = 1> N;               // Number of observations.
  int<lower = 1> K;               // Number of groups.
  matrix[N, K] Y;                 // individual obseravtions
  vector[K] vec_1;
  matrix[K, K] identity_mat;
  matrix[K, K] P_alpha;
  matrix[K, K] P_beta;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> beta;
  real mu_0;
  real sigma;
  vector[K] mu;
  cov_matrix[K] cov_mat; 
}

model {
  mu_0 ~ normal(0, 1);
  sigma ~ inv_gamma(1, 1);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  mu ~ multi_normal(mu_0 * vec_1, sigma * identity_mat);
  cov_mat ~ inv_wishart(100, alpha * P_alpha + beta * P_beta);
  
  Y[K] ~ multi_normal(mu, cov_mat);
}

generated quantities {
  vector[K] y_rep = multi_normal_rng(mu, cov_mat);
}

