set.seed(42)

# ===== Simulated Inputs =====
n_sites <- 100
n_years <- 2
n_visits_np <- 5

# Covariates
z <- rnorm(n_sites)
x <- 0.3 * z + sqrt(1 - 0.3^2) * rnorm(n_sites)  # x ~ correlated with z
x_mat <- matrix(rep(x, n_years), nrow = n_sites)

# Simulated detection histories
D_np <- array(rbinom(n_sites * n_years * n_visits_np, 1, 0.7),
              dim = c(n_sites, n_years, n_visits_np))
D_p <- array(rbinom(n_sites * n_years, 1, 0.7),
             dim = c(n_sites, n_years, 1))

# Observed presence/absence
obs_np <- apply(D_np, c(1,2), max)
obs_p  <- apply(D_p,  c(1,2), max)

# Year matrix: 0 for year 1, 1 for year 2
year_mat <- matrix(rep(0:(n_years - 1), each = n_sites), nrow = n_sites)

# ===== Likelihood Function =====
loglik_fn <- function(par) {
  beta0_p   <- par[1]
  beta0_np  <- par[2]
  beta_year <- par[3]  # shared
  beta_x_p  <- par[4]
  beta_x_np <- par[5]
  
  ll <- 0
  
  # Probability sample
  for (i in 1:n_sites) {
    for (t in 1:n_years) {
      eta_p <- beta0_p + beta_year * year_mat[i,t] + beta_x_p * x[i]
      psi_p <- plogis(eta_p)
      y <- obs_p[i,t]
      ll <- ll + y * log(psi_p) + (1 - y) * log(1 - psi_p)
    }
  }
  
  # Non-probability sample
  for (i in 1:n_sites) {
    for (t in 1:n_years) {
      eta_np <- beta0_np + beta_year * year_mat[i,t] + beta_x_np * x[i]
      psi_np <- plogis(eta_np)
      y <- obs_np[i,t]
      ll <- ll + y * log(psi_np) + (1 - y) * log(1 - psi_np)
    }
  }
  
  return(-ll)  # negative log-likelihood
}

# ===== Optimization =====
start_vals <- rep(0, 5)
fit <- optim(start_vals, loglik_fn, method = "BFGS", hessian = TRUE)

# Output
mle <- fit$par
names(mle) <- c("beta0_p", "beta0_np", "beta_year", "beta_x_p", "beta_x_np")
se <- sqrt(diag(solve(fit$hessian)))

list(mle = mle, se = se)

# Confidence intervals: Wald-type (estimate ± 1.96 * SE)
ci_lower <- mle - 1.96 * se
ci_upper <- mle + 1.96 * se

conf_int <- data.frame(
  Parameter = names(mle),
  Estimate = round(mle, 3),
  SE = round(se, 3),
  CI_Lower = round(ci_lower, 3),
  CI_Upper = round(ci_upper, 3)
)

print(conf_int)

# Empirical marginal occupancy in each year
occ_year0 <- mean(occupancy_matrix[, 1])
occ_year1 <- mean(occupancy_matrix[, 2])

# Clip if needed to avoid 0 or 1 in logit
clip <- function(p) pmax(pmin(p, 0.999), 0.001)

# Log-odds ratio of occupancy (year 1 vs year 0)
logit <- function(p) log(p / (1 - p))

psi0 <- clip(occ_year0)
psi1 <- clip(occ_year1)

true_beta_year <- logit(psi1) - logit(psi0)

cat("True beta_year (from simulated occupancy change):", round(true_beta_year, 3), "\n")

