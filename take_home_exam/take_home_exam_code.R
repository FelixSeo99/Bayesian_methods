
library(mvtnorm)
library(tidyverse)
library(rstan)
library(knitr)


# loading the data
data_toxic <- read_delim("data/toxic.csv", delim = ",")
data_returns <- read_delim("data/returns.csv", delim = ",") %>% 
  mutate(across(contains("Stock"), ~as.numeric(.x)))

# Task 5 b

mean_vec <- c(0, 10)
cov_mat <- matrix(c(4, 12, 12, 100), nrow = 2)

# x and y is vector of values found in the data file. n is the batch size 
# which is 5 here. Outputs the unnormalized log-posterior. 
unnorm_post <- function(alpha, beta, x, y, n){
  prob_pi <- 1 / (1 + exp( - (alpha + beta * x)))
  log_lik <- dbinom(y, n, prob_pi, log = TRUE) # log=TRUE gives the log like.
  log_prior <- dmvnorm(c(alpha, beta), mean_vec, cov_mat, log = TRUE)
  log_post <- sum(log_lik) + log_prior
  log_post 
}

seq_alpha <- seq(-2.5, 5, 0.1)
seq_beta <- seq(-1, 30, 1)


data.frame(
  "alpha" = rep(seq_alpha, each = length(seq_beta)),
  "beta" = rep(seq_beta, times = length(seq_alpha))
) %>% 
  mutate(
    un_post = mapply(function(x, y) unnorm_post(x, y, data_toxic$x, data_toxic$y, 5),
                     alpha, 
                     beta)
  ) %>% # our function not vectorized, not work with mutate.
  ggplot(aes(x = alpha, y = beta, z = un_post)) +
  geom_contour_filled() +
  xlab(latex2exp::TeX(r'($\alpha$)')) +
  ylab(latex2exp::TeX(r'($beta$)')) +
  theme(legend.key.height= unit(0.45, 'cm'), axis.title.y = element_text(angle = 0, vjust = 0.5)) 


# Task 5 c

# Metropolis-Hastings algorithm:
# input starting values of alpha and beta in alpha_init and beta_init 
# respectively. x and y is the data and n the batch size. Iterations is how many
# iterations to go for. Outputs data frame with sequence.
MH_alg <- function(alpha_init, beta_init, x, y, n = 5, iter){
  prop_cov <- matrix(c(1, 0, 0, 5), nrow = 2)
  #theta_seq <- c(theta_init)
  #theta <- theta_init     
  alpha <- alpha_init
  beta <- beta_init
  alpha_seq <- c(alpha_init)
  beta_seq <- c(beta_init)
  
  for (i in 1:iter) {
    prop_sample <- rmvnorm(1, c(alpha, beta), prop_cov) %>% as.vector() # vector easy to work with
    prop_theta <- dmvnorm(c(alpha, beta), prop_sample, prop_cov) # q(theta^{t-1} | theta^{star})
    prop_theta_star <- dmvnorm(prop_sample, c(alpha, beta), prop_cov) # q(theta^{star} | theta^{t-1})
    
    post_theta_star <- unnorm_post(prop_sample[1], prop_sample[2] , x, y, n)
    post_theta <- unnorm_post(alpha, beta, x, y, n) # f(theta^{t-1} | x)
    
    post_fraction <- exp(post_theta_star - post_theta)
    prop_fraction <- prop_theta / prop_theta_star
    accept_prob <- min(1, post_fraction * prop_fraction) # multiplication since calc prop above.
    
    unif <- runif(1, 0, 1)
    if (unif <= accept_prob) {
      alpha <- prop_sample[1]
      beta <- prop_sample[2]
      alpha_seq <- c(alpha_seq, prop_sample[1])
      beta_seq <- c(beta_seq, prop_sample[2])
    } else if (unif > accept_prob) {
      alpha_seq <- c(alpha_seq, alpha)
      beta_seq <- c(beta_seq, beta)
    } 
  }
  
  data.frame("alpha" = alpha_seq, "beta" = beta_seq)
}


# Task 5d traceplots 

set.seed(990108) 

samples <- MH_alg(-2.5, 0, data_toxic$x, data_toxic$y, iter = 10000) 

samples2 <- MH_alg(0, 2.5, data_toxic$x, data_toxic$y, iter = 10000) %>% 
  mutate(group = 2)
samples3 <- MH_alg(1, 7, data_toxic$x, data_toxic$y, iter = 10000) %>% 
  mutate(group = 3)
samples4 <- MH_alg(2.5, 13, data_toxic$x, data_toxic$y, iter = 10000) %>% 
  mutate(group = 4)
samples5 <- MH_alg(4, 20, data_toxic$x, data_toxic$y, iter = 10000) %>% 
  mutate(group = 5)
samples6 <- MH_alg(5, 26, data_toxic$x, data_toxic$y, iter = 10000) %>% 
  mutate(group = 6)

samples %>% 
  mutate(group = 1) %>% 
  bind_rows(samples2, samples3, samples4, samples5, samples6) %>%
  mutate(time = rep(seq(0, 10000), times = 6)) %>% 
  pivot_longer(cols = c("alpha", "beta"), names_to = "parameter", values_to = "value") %>% 
  mutate(group = as.factor(group)) %>% 
  ggplot(aes(x = time, y = value, color = group)) +
  geom_line() +
  facet_wrap(~parameter) +
  labs(color = "Chain")




# Task 5d marginal distribution of the parameters 

samples %>% 
  pivot_longer(cols = everything(), names_to = c("parameter"), values_to = "value") %>% 
  ggplot(aes(x = value, y = after_stat(density))) +
  geom_histogram(bins = 80, color = "black", fill = "white") + 
  facet_wrap(vars(parameter), scales = "free_x")


# task 5d plot sample on the contour plot

test_alpha <- samples3 %>% select(alpha) %>% slice_sample(n = 2432) %>% {.$alpha}
test_beta <- samples3 %>% select(beta) %>% slice_sample(n = 2432) %>% {.$beta}

data.frame(
  "alpha" = rep(seq_alpha, each = length(seq_beta)),
  "beta" = rep(seq_beta, times = length(seq_alpha))
) %>% 
  mutate(
    un_post = mapply(function(x, y) unnorm_post(x, y, data_toxic$x, data_toxic$y, 5),
                     alpha, 
                     beta)
  ) %>% # our function not vectorized, not work with mutate.
  ggplot(aes(x = alpha, y = beta, z = un_post)) +
  geom_contour_filled() +
  geom_point(data = NULL, aes(x = test_alpha, y = test_beta)) +
  xlab(latex2exp::TeX(r'($\alpha$)')) +
  ylab(latex2exp::TeX(r'($beta$)')) +
  theme(legend.key.height= unit(0.45, 'cm'), axis.title.y = element_text(angle = 0, vjust = 0.45)) 

########################


# Task 6


# Task 6a the summary statistics for the data

data_returns_2 <- data_returns %>% 
  select(-Date) %>% 
  pivot_longer(everything(), names_to = "stock", values_to = "returns") 

data_returns_2 %>% 
  group_by(stock) %>% 
  summarise(
    mean = mean(returns), 
    variance = var(returns), 
    median = median(returns)
  ) %>% 
  kable(
    align = "c", 
    caption = "Showing the mean, variance and median of the returns of the 6 different stocks.")


# Boxplots of the data

data_returns_2 %>% 
  ggplot(aes(x = stock, y = returns)) +
  geom_boxplot()

# Histograms for the stocks  

data_returns_2 %>% 
  ggplot(aes(x = returns, y = after_stat(density))) +
  geom_histogram(bins = 50) + 
  facet_wrap(~ stock)


# Task 6 b

# Define the matrices P_A and P_B. These are matrices where we only have 
# the P_A block in the 6x6 matrix P_0 and not multiplied by alpha. all other 
# elements are 0. Same for P_B but for the block P_B. See stan file for more 
# info on how we then get the P_0 matrix. 
P_A <- matrix(c(1, 0.5, 0.5, 0, 0, 0, 
                0.5, 1, 0.5, 0, 0, 0, 
                0.5, 0.5, 1, 0, 0, 0,
                rep(0, 18)
),
nrow = 6) 

P_B <- matrix(c(0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0.5, 0.5, 
                0, 0, 0, 0.5, 1, 0.5,
                0, 0, 0, 0.5, 0.5, 1
),
nrow = 6)

data_stock <- data_returns %>% select(!Date) %>% as.matrix()

model_data <- list(
  "N" = length(data_returns$Date),
  "K" = 6,
  "Y" = data_stock,
  "vec_1" = rep(1, 6),
  "identity_mat" = diag(1, 6), 
  "P_A" = P_A, 
  "P_B" = P_B
)

model_fit <- stan(
  "take_home_exam/stan_model_6.stan",
  data = model_data,
  chains = 4, 
  iter = 30000, 
  warmup = 20000,
  seed = 990108,
  cores = 4
)

# Task 6d

# Deriving summary statistics for the predictive checks. 
rstan::extract(model_fit, pars = "y_rep") %>% 
  as.data.frame() %>% 
  rename(
    "Energy.Stock1" = "y_rep.1",
    "Energy.Stock2" = "y_rep.2",
    "Energy.Stock3" = "y_rep.3",
    "IT.Stock1" = "y_rep.4",
    "IT.Stock2" = "y_rep.5",
    "IT.Stock3" = "y_rep.6"
  ) %>% 
  pivot_longer(everything(), names_to = "stock", values_to = "returns") %>% 
  group_by(stock) %>% 
  summarise(mean = mean(returns), variance = var(returns), median = median(returns)) %>% 
  kable(
    align = "c", 
    caption = "Showing the mean, variance and median based on a sample from the 
    posterior predictive distribution.")


# Traceplots given in the same order as they appear in the report.

traceplot(model_fit, pars =c("alpha", "beta", "sigma_2", "mu_0")) + 
  xlab("iteration")

traceplot(model_fit, pars = "mu") + 
  xlab("iteration")


# the rhat values for the mean vector elements. 
summary(model_fit)$summary[,"Rhat"] %>% 
  as.data.frame() %>% 
  rownames_to_column("parameter") %>% 
  rename("Rhat" = ".") %>% 
  filter(grepl("mu", parameter) & parameter != "mu_0") %>%
  pivot_wider(names_from = "parameter", values_from = "Rhat") %>% 
  kable(
    align = "c",
    col.names = c("$\\mu_1$", "$\\mu_2$", "$\\mu_3$", "$\\mu_4$", "$\\mu_5$", "$\\mu_6$"),
    caption = "Rhat values (Gelman-Rubin statistic) for the elements in the mean vector $\\boldsymbol{\\mu}$.")


# autocorrealtion plot for mean vector
stan_ac(model_fit, pars = "mu")


# traceplot diagonal elements in covariance matrix
traceplot(
  model_fit, 
  pars = c(
    "cov_mat[1,1]", 
    "cov_mat[2,2]", 
    "cov_mat[3,3]",
    "cov_mat[4,4]", 
    "cov_mat[5,5]",
    "cov_mat[6,6]"
  )
) + 
  xlab("iteration")


# autocorrelation plot for diagonal elements in covariance matrix
stan_ac(
  model_fit, 
  pars = c(
    "cov_mat[1,1]", 
    "cov_mat[2,2]", 
    "cov_mat[3,3]",
    "cov_mat[4,4]", 
    "cov_mat[5,5]",
    "cov_mat[6,6]"
  )
) 

# Traceplots for the upper non-diagonal elements corresponding to covariances
# between the stocks, i.e. the upper non-diagonal elements in the covariance
# matrix $\\boldsymbol{\\Sigma}$.
traceplot(
  model_fit, 
  pars = c(
    "cov_mat[1,4]", 
    "cov_mat[1,5]", 
    "cov_mat[1,6]",
    "cov_mat[2,4]", 
    "cov_mat[2,5]",
    "cov_mat[2,6]",
    "cov_mat[3,4]", 
    "cov_mat[3,5]",
    "cov_mat[3,6]"
  )
) + 
  xlab("iteration")

# Traceplots of the non-diagonal elements in the covariance matrix 
# $\\boldsymbol{\\Sigma}$ corresponding to the non-diagonal elements in the 
# $\\alpha\\textbf{P}_A$ block of the matrix $\\textbf{P}_0$.

traceplot(
  model_fit, 
  pars = c(
    "cov_mat[1,2]", 
    "cov_mat[1,3]", 
    "cov_mat[2,3]"
  )
) + 
  xlab("iteration")



