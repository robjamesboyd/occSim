find_sampling_intercept <- function(target_f, Z, beta, V = 1, tol = 1e-5, maxit = 500) {
  # f(α) = expected sampling fraction for a given α
  f_alpha <- function(alpha) {
    delta <- plogis(alpha + beta * Z)
    mean(1 - (1 - delta)^V)
  }
  
  # Define the root function: difference between actual and target
  root_func <- function(alpha) f_alpha(alpha) - target_f
  
  # Use uniroot to solve
  out <- uniroot(root_func, interval = c(-20, 20), tol = tol, maxiter = maxit)
  
  return(out$root)
}

target_f <- 0.1
Z <- z
V <- 1
beta <- 2
alpha <- find_sampling_intercept(target_f=0.1, Z=z, V=1,beta=2)
alpha
