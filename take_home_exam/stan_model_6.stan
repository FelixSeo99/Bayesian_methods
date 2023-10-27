// Insert thh number of observation N, number of stocks K, observations Y,
// an identity matrix size 6, a vector of 1's of length 6, and P_A and 
// P_B. The P_B corresponds to the P_A block in the P_0 matrix but 
// all other blocks are 0. So we have a 6x6 matrix with only the P_A block and
// not multiplied by alpha. The same for P_B but here we have the P_B block only
// in a 6x6 matrix. The multiplication with alpha and beta is accounted for 
// later in the below code and also that we get the correct amtrix P_0. 


data {
  int<lower = 1> N;               // Number of observations.
  int<lower = 1> K;               // Number of stocks.
  matrix[N, K] Y;                 // Individual obseravtions
  vector[K] vec_1;                // vector of 1's
  matrix[K, K] identity_mat;      // identity matrix size 6  
  matrix[K, K] P_A;               
  matrix[K, K] P_B;
}

parameters {
  real<lower = 0> alpha;          // set to accept non-negative values
  real<lower = 0> beta;
  real mu_0;
  real sigma_2;                   // sigma^2      
  vector[K] mu;
  cov_matrix[K] cov_mat; 
}

model {
  mu_0 ~ normal(0, 1);
  sigma_2 ~ inv_gamma(1, 1);      // sigma^2
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  mu ~ multi_normal(mu_0 * vec_1, sigma_2 * identity_mat);
  cov_mat ~ inv_wishart(100, alpha * P_A * 100 + beta * P_B * 100);  // sum to 
                                                                     // get P_0
  Y[K] ~ multi_normal(mu, cov_mat);
}

// Generates a sample from the posterior predictive distribution
generated quantities {
  vector[K] y_rep;
  y_rep = multi_normal_rng(mu, cov_mat);
}

