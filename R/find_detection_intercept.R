find_detection_intercept <- function(target_error, w, beta, V = 1, tol = 1e-6, maxit = 100) {
  f_alpha <- function(alpha) {
    tau <- plogis(alpha + beta * w)
    mean((1 - tau)^V)
  }
  
  root_func <- function(alpha) f_alpha(alpha) - target_error
  
  out <- uniroot(root_func, interval = c(-20, 20), tol = tol, maxiter = maxit)
  return(out$root)
}
