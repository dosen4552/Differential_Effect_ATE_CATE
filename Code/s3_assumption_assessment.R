library(foreign)
library(dplyr)
library(purrr)
library(rstan)



#### Implementation of real data ####

################ Lead ##########################
NHANES <- read.csv("kan_bingkai_dylan/Differential_Effect/Data/NHANES/NHANES_2PL_new.csv")
data_example_lead_age <- NHANES
colnames(data_example_lead_age)[1:3] <- c('Y',"Treatment",'Z')


############### 2PL model ################
# convert race to factor
data_example_lead_age$race <- as.factor(data_example_lead_age$race)
data_example_lead_age$education <- as.factor(data_example_lead_age$education)

# Convert categorical variable  to one-hot encoding
race_onehot <- model.matrix(~ race - 1, data = data_example_lead_age)
race_onehot_numeric <- as.data.frame(data.matrix(race_onehot))

# educ_onehot <- model.matrix(~ education - 1, data = data_example_lead_age)
# educ_onehot_numeric <- as.data.frame(data.matrix(educ_onehot))

covariates <- data_example_lead_age |> dplyr::select(age, gender, education,
                                    income_ratio) |> cbind(race_onehot_numeric)

items <- data_example_lead_age %>% dplyr::select(Treatment, Z, fat_intake, alcohol_abuse)
items <- as.data.frame(lapply(items, as.numeric))

# setup stan stuff
N <- nrow(data_example_lead_age)  # Number of respondents
J <- ncol(items)    # Number of items
P <- ncol(covariates)    # Number of covariates
X <- covariates
Y <- items

stan_data <- list(
  N = N,
  J = J,
  P = P,
  X = X,
  Y = Y
)

# Compile the Stan model
stan_model <- stan_model("kan_bingkai_dylan/Differential_Effect/Code/binary_latent_nhanes.stan")

# Fit the model
#fit <- sampling(
#  stan_model,
#  data = stan_data,
#  chains = 4,
#  cores = 4,
#  iter = 2000,
#  warmup = 1000,
#  seed = 123
#)

# Fit the model
fit <- sampling(
  stan_model,
  data = stan_data,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

posterior <- extract(fit)

# Alpha estimates
alpha_samples <- posterior$alpha
alpha_means <- colMeans(alpha_samples)
print(alpha_means)  # Posterior means for alpha
alpha_medians <- apply(alpha_samples, 2, median)
print(alpha_medians)

alpha_sd <- apply(alpha_samples, 2, sd)
print(alpha_sd)


# Beta estimates
beta_samples <- posterior$beta
beta_means <- apply(beta_samples, c(2,3), mean)
print(beta_means)  # Posterior means for alpha
beta_medians <- apply(beta_samples, c(2,3), median)
print(beta_medians)


# Beta estimates
beta_0_samples <- posterior$beta_0
beta_0_means <- colMeans(beta_0_samples)
print(beta_0_means)
beta_0_medians <- apply(beta_0_samples, 2, median)
print(beta_0_medians)


# use posterior means as point estimates
alpha_1 <- alpha_means[1]
alpha_2 <- alpha_means[2]
beta_1 <- beta_0_means[1]
beta_2 <- beta_0_means[2]
cov_1 <- beta_means[,1] 
cov_2 <- beta_means[,2]

# use posterior medians as point estimates
alpha_1 <- alpha_medians[1]
alpha_2 <- alpha_medians[2]
beta_1 <- beta_0_medians[1]
beta_2 <- beta_0_medians[2]
cov_1 <- beta_medians[,1] 
cov_2 <- beta_medians[,2] 


# Initialize an empty list to store results
results_10 <- rep(NA, N)
results_01 <- rep(NA, N)

# Loop through each row in covariates
for (i in 1:nrow(covariates)) {
  current_row <- as.numeric(covariates[i, ]) # Extract the current row as a numeric vector
  
  # Apply the function for all values of u
  results_for_row_10 <- check_prob_diff_10(
        alpha_1, alpha_2, beta_1, beta_2, cov_1, cov_2, current_row
      )

  results_for_row_01 <- check_prob_diff_01(
    alpha_1, alpha_2, beta_1, beta_2, cov_1, cov_2, current_row
  )
  
  # Append the results to the list
  results_10[i] <- results_for_row_10
  results_01[i] <- results_for_row_01
}

hist(results_01, breaks = 1000)
hist(results_10, breaks = 1000)

# proportion of positive sign 
mean(results_01 > 0)
# proportion of negative sign
mean(results_10 < 0)


# functions to compute probabilities and their differences
check_prob_01 <- function(alpha_1, alpha_2, beta_1,
                          beta_2, cov_1, cov_2, u, covariates) {
  exp(alpha_2 * u + beta_2 +sum(cov_2 * covariates)) /
    (1 +  exp(alpha_2 * u + beta_2 +sum(cov_2 * covariates)) + 
       exp(alpha_1 * u + beta_1 +sum(cov_1 * covariates)) + 
       exp(alpha_1 * u + beta_1 +sum(cov_1 * covariates) + 
             alpha_2 * u + beta_2 +sum(cov_2 * covariates)))
}

check_prob_diff_01 <- function(alpha_1, alpha_2, beta_1,
                               beta_2, cov_1, cov_2, covariates) {
  check_prob_01(alpha_1, alpha_2, beta_1,
                beta_2, cov_1, cov_2, 1, covariates) - 
    check_prob_01(alpha_1, alpha_2, beta_1,
                  beta_2, cov_1, cov_2, 0, covariates)
}

check_prob_10 <- function(alpha_1, alpha_2, beta_1,
                          beta_2, cov_1, cov_2, u, covariates) {
  exp(alpha_1 * u + beta_1 +sum(cov_1 * covariates)) /
    (1 +  exp(alpha_2 * u + beta_2 +sum(cov_2 * covariates)) + 
       exp(alpha_1 * u + beta_1 +sum(cov_1 * covariates)) + 
       exp(alpha_1 * u + beta_1 +sum(cov_1 * covariates) + 
             alpha_2 * u + beta_2 +sum(cov_2 * covariates)))
}

check_prob_diff_10 <- function(alpha_1, alpha_2, beta_1,
                               beta_2, cov_1, cov_2, covariates) {
  check_prob_10(alpha_1, alpha_2, beta_1,
                beta_2, cov_1, cov_2, 1, covariates) - 
    check_prob_10(alpha_1, alpha_2, beta_1,
                  beta_2, cov_1, cov_2, 0, covariates)
}


