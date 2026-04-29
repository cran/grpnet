row.kronecker <- 
  function(X, Y){
    ### row-wise kronecker product
    ### Nathaniel E. Helwig (helwig@umn.edu)
    ### Mar 31, 2026
    
    # check dimensions
    nx <- nrow(X)
    ny <- nrow(Y)
    if(nx != ny) stop("Inputs 'X' and 'Y' must have the same number of rows")
    px <- ncol(X)
    py <- ncol(Y)
    pz <- px * py
    
    # check names
    xnames <- colnames(X)
    if(is.null(xnames)) xnames <- paste0("X", 1:px)
    ynames <- colnames(Y)
    if(is.null(ynames)) ynames <- paste0("Y", 1:py)
    
    # initialize and build Z
    znames <- vector(mode = "character", length = pz)
    Z <- matrix(0.0, nrow = nx, ncol = pz)
    for(j in 1:px){
      index <- 1:py + (j - 1) * py
      Z[,index] <- X[,j] * Y
      znames[index] <- paste0(xnames[j], ":", ynames)
    }
    colnames(Z) <- znames
    return(Z)
    
  }