R_grpnet_maxeigval <- 
  function(A, N, MAXEV){
    A <- as.matrix(A)
    N <- as.integer(N)
    X <- rep(1 / sqrt(N), N)
    MAXEV <- 1.0
    OLDEV <- 0.0
    while(abs(MAXEV - OLDEV) >= 1e-8){
      OLDEV <- MAXEV
      X <- A %*% X
      MAXEV <- sqrt(sum(X^2))
      if(MAXEV > 0.0) X <- X / MAXEV
    }
    return(MAXEV)
  } # grpnet_maxeigval