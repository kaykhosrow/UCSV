UCSV <- function(x, T, nloop, burnin, Vtau, Vh, atau, ltau, ah, lh) {
  # Initialize the Markov chain
  omega2tau = .25^2
  omega2h = .2^2
  h = as.numeric(log(var(y) * .8)) * matrix(rep(1, T))
  
  # Creating the sparse matrix
  Sp = diag(T)
  d = matrix(rep(0, T), 1)
  sparse = rbind(d, Sp)
  sparse = sparse[-(T + 1), ]
  H = diag(T) - sparse
  
  # Initialize storage
  store_omega2tau = matrix(0, nloop - burnin, 1)
  store_omega2h = matrix(0, nloop - burnin, 1)
  store_tau = matrix(0, (nloop - burnin) * T, T)
  store_h = matrix(0, (nloop - burnin) * T, T)
  
  # Precompute some constants
  newatau = (T - 1) / 2 + atau
  newah = (T - 1) / 2 + ah
  
  for (loop in 1:nloop) {
    # Sample tau
    invOmegatau = diag(T) * c(1 / Vtau, 1 / omega2tau * rep(1, T - 1))
    invSigy = diag(T) * exp(-h)
    Ktau = t(H) %*% invOmegatau %*% H + invSigy
    Ctau = t(chol(Ktau))
    tauhat = solve(Ktau, invSigy %*% y)
    tau = tauhat + solve(t(Ctau), matrix(rnorm(T), 1))
    
    # Sample h
    ystar = log((y - tau)^2 + .0001)
    result = support(ystar, h, omega2h, Vh)
    h = result[[1]]
    
    # Sample omega2tau
    newltau = ltau + sum((tau[2:T] - tau[1:(T - 1)])^2) / 2
    omega2tau = gamrand(alpha = newatau, lambda = newltau)
    
    # Sample omega2h
    newlh = lh + sum((h[2:T] - h[1:(T - 1)])^2) / 2
    omega2h = gamrand(alpha = newah, lambda = newlh)
    
    # Store the results
    if (loop > burnin) {
      i = loop - burnin
      store_tau[i, ] = t(tau)
      store_h[i, ] = t(h)
      store_omega2tau[i, ] = omega2tau
      store_omega2h[i, ] = omega2h
    }
  }
  
  tauhat = rowMeans(store_tau)
  hhat = colMeans(store_h)
  
  return(list(tauhat = tauhat, hhat = hhat, store_tau = store_tau, store_h = store_h))
}

# Gamma random variable generator function
gamrand <- function(alpha, lambda) {
  if (alpha > 1) {
    d = alpha - 1/3
    c = 1 / sqrt(9 * d)
    flag = TRUE
    while (flag) {
      Z = rnorm(1)
      if (Z > (-1 / c)) {
        V = (1 + c * Z)^3
        U = runif(1)
        flag = log(U) > (0.5 * Z^2 + d - d * V + d * log(V))
      }
    }
    x = d * V / lambda
  } else {
    x = gamrand(alpha + 1, lambda) * runif(1)^(1 / alpha)
  }
  return(x)
}

# Example usage with parameterisation
y = x
T = length(x)
nloop = 11000
burnin = 1000
Vtau = 9
Vh = 9
atau = 10
ltau = .25^2 * (atau - 1)
ah = 10
lh = .2^2 * (ah - 1)

result = UCSV(y, T, nloop, burnin, Vtau, Vh, atau, ltau, ah, lh)

UCSV <- result$tauhat