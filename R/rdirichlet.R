# Simulate values from a Dirichlet distribution
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
  x <- x / rowSums(x)
  if (n == 1) {
    return(as.vector(x))
  } else {
    return(x)
  }
}
