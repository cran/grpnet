R_grpnet_ordinal <-
  function(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw, 
           alpha, nlam, lambda, lmr, penid, gamma, eps, maxit,
           standardize, intercept, ibeta, betas, iters,
           nzgrps, nzcoef, edfs, devs, nulldev){
    # grpnet_ordinal.f90 translation to R
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2026-04-09
    
    
    # ! --------------- LOCAL DEFINITIONS --------------- ! #
    ia <- ib <- rep(0L, ngrps)
    xmean <- rep(0.0, nvars)
    xev <- rep(0.0, ngrps)
    difbeta <- rep(0.0, nvars)
    # ! --------------- LOCAL DEFINITIONS --------------- ! #
    
    
    # ! --------------- MACHINE EPSILON --------------- ! #
    macheps <- .Machine$double.eps
    # ! --------------- MACHINE EPSILON --------------- ! #
    
    
    # ! --------------- CHECK WEIGHTS --------------- ! #
    wmin <- min(w)
    wmax <- max(w)
    if(wmax > wmin){
      weighted <- 1L
      w <- nobs * w / sum(w)     # ! normalize so SUM(w) = nobs
      w <- sqrt(w)
      for(i in 1:nobs){
        x[i,] <- w[i] * x[i,]
      }
    } else {
      weighted <- 0L
      w <- rep(1.0, nobs)
    }
    # ! --------------- CHECK WEIGHTS --------------- ! #
    
    
    # ! --------------- GROUP INDICES --------------- ! #
    gid <- 0L
    for(k in 1:ngrps){
      ia[k] <- gid + 1L
      gid <- gid + gsize[k]
      ib[k] <- gid
    }
    # ! --------------- GROUP INDICES --------------- ! #
    
    
    # ! --------------- CENTER AND SCALE --------------- ! #
    if(intercept == 1L){
      for(j in 1:nvars){
        xmean[j] <- sum(x[,j]) / nobs
        x[,j] <- x[,j] - xmean[j]
      }
    }
    xsdev <- rep(1.0, ngrps)
    if(standardize == 1L){
      for(k in 1:ngrps){
        xsdev[k] <- sqrt( sum(x[,ia[k]:ib[k]]^2) / (nobs * gsize[k]) )
        if(xsdev[k] > macheps){
          x[,ia[k]:ib[k]] <- x[,ia[k]:ib[k]] / xsdev[k]
        } else {
          xsdev[k] <- 1.0
        }
      }
    }
    # ! --------------- CENTER AND SCALE --------------- ! #
    
    
    # ! --------------- GET MAX EIGENVALUE --------------- ! #
    for(k in 1:ngrps){
      if(gsize[k] == 1L){
        xev[k] <- sum(x[,ia[k]]^2) / nobs
      } else {
        xtx <- crossprod(x[,ia[k]:ib[k]]) / nobs
        xev[k] <- R_grpnet_maxeigval(A = xtx, N = gsize[k])
      }
    }
    xev <- xev / 2.0
    # ! --------------- GET MAX EIGENVALUE --------------- ! #
    
    
    # ! --------------- MISCELLANEOUS INITIALIZATIONS --------------- ! #
    maxlam <- max(lambda)
    makelambda <- 0L
    iter <- 0L
    active <- strong <- rep(0L, ngrps)
    nzgrps <- nzcoef <- rep(0L, nlam)
    ibeta <- matrix(0.0, nresp-1L, nlam)
    beta <- zvec <- grad <- rep(0.0, nvars)
    gradnorm <- rep(0.0, ngrps)
    twolam <- 0.0
    devs <- rep(0.0, nlam)
    eta <- off
    # ! --------------- MISCELLANEOUS INITIALIZATIONS --------------- ! #
    
    
    # ! --------------- SHIFT INITIALIZATION --------------- ! #
    yfreq <- ycumsum <- rep(0L, nresp)
    yfreq[1] <- ycumsum[1] <- sum(y == 1L)
    for(k in 2:nresp) {
      yfreq[k] <- sum(y == k)
      ycumsum[k] <- ycumsum[k-1] + yfreq[k]
    }
    yprob <- 1 - (ycumsum / nobs)[-nresp]
    yprob[yprob < 1e-6] <- 1e-6
    yprob[yprob > 1 - 1e-6] <- 1 - 1e-6
    shift <- log(yprob / (1 - yprob))
    # ! --------------- SHIFT INITIALIZATION --------------- ! #
    
    
    # ! --------------- GENERATE LAMBDA --------------- ! #
    if(maxlam <= macheps){
      
      makelambda <- 1L
      i <- 1L
      
      # ! find unpenalized groups ! # 
      for(k in 1:ngrps){
        if(pw[k] <= macheps){
          active[k] <- 1L
          nzgrps[i] <- nzgrps[i] + 1L
          nzcoef[i] <- nzcoef[i] + gsize[k]
        }
      }
      # ! find unpenalized groups ! #
      
      # ! iterate until active coefficients converge ! #
      if(nzgrps[i] > 0L){
        
        while(iter < maxit){
          
          # ! update iter and reset counters
          ctol <- 0.0
          iter <- iter + 1L
          
          # ! update shifts
          temp <- R_grpnet_ordinal_int(nobs, nresp, eta, y, w^2, yfreq, shift, 1.0)
          shift <- temp$theta
          ctol <- max(temp$maxdif, ctol)
          
          # ! define r
          bigtheta <- c(100, shift, -100)
          fit0 <- eta + bigtheta[y]
          fit1 <- eta + bigtheta[y+1L]
          Pi0 <- 1 / (1 + exp(-fit0))
          Pi1 <- 1 / (1 + exp(-fit1))
          r <- w * (1.0 - Pi0 - Pi1)
          
          # ! update active groups
          for(k in 1:ngrps){
            if(active[k] == 0L) next
            grad[ia[k]:ib[k]] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
            zvec[ia[k]:ib[k]] <- beta[ia[k]:ib[k]] + grad[ia[k]:ib[k]] / xev[k]
            difbeta[ia[k]:ib[k]] <- zvec[ia[k]:ib[k]] - beta[ia[k]:ib[k]]
            maxdif <- max( abs(difbeta[ia[k]:ib[k]]) / (1.0 + abs(beta[ia[k]:ib[k]])) )
            beta[ia[k]:ib[k]] <- beta[ia[k]:ib[k]] + difbeta[ia[k]:ib[k]]
            eta <- eta + (x[,ia[k]:ib[k]] %*% difbeta[ia[k]:ib[k]]) / w
            fit0 <- eta + bigtheta[y]
            fit1 <- eta + bigtheta[y+1L]
            Pi0 <- 1 / (1 + exp(-fit0))
            Pi1 <- 1 / (1 + exp(-fit1))
            r <- w * (1.0 - Pi0 - Pi1)
            ctol <- max(maxdif, ctol)
          } # ! k=1,ngrps
          
          # ! convergence check
          if(ctol < eps) break
          
        } # ! WHILE(iter < maxit)
        
      } else {
        
        # ! intercept only
        while(iter < maxit){
          
          # ! update iter and reset counters
          ctol <- 0.0
          iter <- iter + 1L
          
          # ! update shifts
          temp <- R_grpnet_ordinal_int(nobs, nresp, eta, y, w^2, yfreq, shift, 1.0)
          shift <- temp$theta
          
          # ! convergence check
          ctol <- max(temp$maxdif, ctol)
          if(ctol < eps) break
          
        }  # ! WHILE(iter < maxit)
        
        # ! define r
        bigtheta <- c(100, shift, -100)
        fit0 <- eta + bigtheta[y]
        fit1 <- eta + bigtheta[y+1L]
        Pi0 <- 1 / (1 + exp(-fit0))
        Pi1 <- 1 / (1 + exp(-fit1))
        r <- w * (1.0 - Pi0 - Pi1)
        
      } # ! (nzgrps[i] > 0)
      # ! iterate until active coefficients converge ! #
      
      # ! create lambda sequence ! #
      for(k in 1:ngrps){
        if(pw[k] > macheps){
          grad[ia[k]:ib[k]] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
          gradnorm[k] <- sqrt( sum(grad[ia[k]:ib[k]]^2) ) / pw[k]
        }
      }
      if(alpha > macheps){
        maxlam <- max(gradnorm / alpha)
      } else {
        maxlam <- max(gradnorm / 1e-3)
        makelambda <- 0L
      }
      minlam <- lmr * maxlam
      lambda[i] <- maxlam
      maxlam <- log(maxlam)
      minlam <- log(minlam)
      rnglam <- maxlam - minlam
      for(k in 2:nlam){
        lambda[k] <- exp( maxlam - rnglam * (k - 1) / (nlam - 1) )
      }
      # ! create lambda sequence ! #
      
      # ! calculate deviance ! #
      mu <- Pi0 - Pi1
      devs[i] <- -2.0 * sum(w^2 * log(mu))
      # ! calculate deviance ! #
      
      # ! save results ! # 
      ibeta[,i] <- shift
      betas[,i] <- beta
      iters[i] <- iter
      nzgrps[i] <- nzgrps[i] + intercept
      nzcoef[i] <- nzcoef[i] + intercept*(nresp-1)
      edfs[i] <- as.numeric(nzcoef[i])
      # ! save results ! # 
      
    }
    # ! --------------- GENERATE LAMBDA --------------- ! #
    
    
    # ! --------------- ITERATIVE WORK --------------- ! #
    for(i in 1:nlam){
      
      # ! initializations ! # 
      if(i == 1L && makelambda == 1L) next
      if(i > 1L){
        shift <- ibeta[,i-1]
        beta <- betas[,i-1]
        twolam <- alpha * (2.0 * lambda[i] - lambda[i-1])
      } else {
        grad <- crossprod(x, r) / nobs
      }
      # ! initializations ! # 
      
      # ! strong rule initialization ! #
      for(k in 1:ngrps){
        gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k]]^2))
        if(gradnorm[k] + 1e-8 > pw[k] * twolam){
          strong[k] <- 1L
        } else {
          strong[k] <- 0L
        }
      }
      # ! strong rule initialization ! #
      
      # ! iterate until strong set converges ! #
      iter <- 0L
      while(iter < maxit){
        
        # ! iterate until active set converges ! # 
        while(iter < maxit){
          
          # ! iterate until active coefficients converge ! #
          while(iter < maxit){
            
            # ! update iter and reset counters
            iter <- iter + 1L
            nzgrps[i] <- 0L
            edfs[i] <- 0.0
            ctol <- 0.0
            
            # ! update shifts
            temp <- R_grpnet_ordinal_int(nobs, nresp, eta, y, w^2, yfreq, shift, 1.0)
            shift <- temp$theta
            ctol <- max(temp$maxdif, ctol)
            
            # ! define r
            bigtheta <- c(100, shift, -100)
            fit0 <- eta + bigtheta[y]
            fit1 <- eta + bigtheta[y+1L]
            Pi0 <- 1 / (1 + exp(-fit0))
            Pi1 <- 1 / (1 + exp(-fit1))
            r <- w * (1.0 - Pi0 - Pi1)
            
            # ! update active groups
            for(k in 1:ngrps){
              if(active[k] == 0L) next
              penone <- alpha * lambda[i] * pw[k] / xev[k]
              pentwo <- (1 - alpha) * lambda[i] * pw[k] / xev[k]
              grad[ia[k]:ib[k]] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
              zvec[ia[k]:ib[k]] <- beta[ia[k]:ib[k]] + grad[ia[k]:ib[k]] / xev[k]
              znorm <- sqrt(sum(zvec[ia[k]:ib[k]]^2))
              bnorm <- sqrt(sum(beta[ia[k]:ib[k]]^2))
              shrink <- R_grpnet_penalty(znorm, penid, penone, pentwo, gamma)
              if(shrink == 0.0 && bnorm == 0.0) next
              difbeta[ia[k]:ib[k]] <- shrink * zvec[ia[k]:ib[k]] - beta[ia[k]:ib[k]]
              maxdif <- max( abs(difbeta[ia[k]:ib[k]]) / (1.0 + abs(beta[ia[k]:ib[k]])) )
              beta[ia[k]:ib[k]] <- beta[ia[k]:ib[k]] + difbeta[ia[k]:ib[k]]
              eta <- eta + (x[,ia[k]:ib[k]] %*% difbeta[ia[k]:ib[k]]) / w
              fit0 <- eta + bigtheta[y]
              fit1 <- eta + bigtheta[y+1L]
              Pi0 <- 1 / (1 + exp(-fit0))
              Pi1 <- 1 / (1 + exp(-fit1))
              r <- w * (1.0 - Pi0 - Pi1)
              ctol <- max(maxdif , ctol)
              if(shrink > 0.0){
                nzgrps[i] <- nzgrps[i] + 1L
                edfs[i] <- edfs[i] + gsize[k] * shrink
              }
            } # ! k=1,ngrps
            
            # ! convergence check
            if(ctol < eps) break
            
          } # ! WHILE(iter < maxit) - inner
          # ! iterate until active coefficients converge ! #
          
          # ! check inactive groups in strong set ! #
          violations <- 0L
          for(k in 1:ngrps){
            if(strong[k] == 0L | active[k] == 1L) next
            grad[ia[k]:ib[k]] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
            gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k]]^2))
            if(gradnorm[k] > alpha * lambda[i] * pw[k]){
              active[k] <- 1L
              violations <- violations + 1L
            }
          } # ! k=1,ngrps
          if(violations == 0L) break
          # ! check inactive groups in strong set ! #
          
        } # ! WHILE(iter < maxit) - middle
        # ! iterate until active set converges ! # 
        
        # ! check groups in weak set ! #
        violations <- 0L
        for(k in 1:ngrps){
          if(strong[k] == 1L) next
          grad[ia[k]:ib[k]] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
          gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k]]^2))
          if(gradnorm[k] + 1e-8 > alpha * lambda[i] * pw[k]){
            strong[k] <- 1
            active[k] <- 1
            violations <- violations + 1L
          }
        }
        if(violations == 0L) break
        # ! check groups in weak set ! #
        
      } # ! WHILE(iter < maxit) - outer
      # ! iterate until strong set converges ! #
      
      # ! calculate nzcoef ! #
      for(k in 1:ngrps){
        if(active[k] == 0L) next
        for(l in 1:gsize[k]){
          if(abs(beta[ia[k] + l - 1]) > macheps){
            nzcoef[i] <- nzcoef[i] + 1L
          }
        }
      }
      # ! calculate nzcoef ! #
      
      # ! calculate deviance ! #
      mu <- Pi0 - Pi1
      devs[i] <- -sum(2 * w^2 * log(mu))
      # ! calculate deviance ! #
      
      # ! save results ! #
      ibeta[,i] <- shift
      betas[,i] <- beta
      iters[i] <- iter
      nzgrps[i] <- nzgrps[i] + intercept
      nzcoef[i] <- nzcoef[i] + intercept
      edfs[i] <- edfs[i] + as.numeric(intercept * (nresp - 1))
      # ! save results ! #
      
    }
    # ! --------------- ITERATIVE WORK --------------- ! #
    
    
    # ! --------------- POST PROCESSING --------------- ! #
    if(standardize == 1L){
      for(k in 1:ngrps){
        betas[ia[k]:ib[k],] <- betas[ia[k]:ib[k],] / xsdev[k]
      }
    }
    if(intercept == 1L){
      xmb <- as.numeric(xmean %*% betas)
      for(i in 1:nlam){
        ibeta[,i] <- ibeta[,i] - xmb[i]
      }
    } 
    if(makelambda == 1L && nzgrps[1] == 0L){
      bigtheta <- c(100, ibeta[,1], -100)
      fit0 <- bigtheta[y]
      fit1 <- bigtheta[y+1L]
      Pi0 <- 1 / (1 + exp(-fit0))
      Pi1 <- 1 / (1 + exp(-fit1))
    } else {
      shift <- log(yprob / (1 - yprob))
      eta <- rep(0.0, nobs)
      while(iter < maxit){
        ctol <- 0.0
        iter <- iter + 1L
        temp <- R_grpnet_ordinal_int(nobs, nresp, eta, y, w^2, yfreq, shift, 1.0)
        shift <- temp$theta
        ctol <- max(temp$maxdif, ctol)
        if(ctol < eps) break
      }  # ! WHILE(iter < maxit)
      bigtheta <- c(100, shift, -100)
      fit0 <- bigtheta[y]
      fit1 <- bigtheta[y+1L]
      Pi0 <- 1 / (1 + exp(-fit0))
      Pi1 <- 1 / (1 + exp(-fit1))
    }
    mu <- Pi0 - Pi1
    nulldev <- -sum(2 * w^2 * log(mu))
    names(xsdev) <- names(pw)
    pw <- xsdev
    # ! --------------- POST PROCESSING --------------- ! #
    
    
    # ! --------------- RETURN RESULTS --------------- ! #
    res <- list(nobs = nobs, 
                nvars = nvars, 
                x = x, 
                y = y, 
                w = w, 
                off = off, 
                ngrps = ngrps, 
                gsize = gsize, 
                pw = pw, 
                alpha = alpha, 
                nlam = nlam, 
                lambda = lambda, 
                lmr = lmr, 
                penid = penid, 
                gamma = gamma, 
                eps = eps, 
                maxit = maxit,
                standardize = standardize, 
                intercept = intercept, 
                ibeta = ibeta, 
                betas = betas, 
                iters = iters,
                nzgrps = nzgrps, 
                nzcoef = nzcoef, 
                edfs = edfs, 
                devs = devs, 
                nulldev = nulldev)
    return(res)
    # ! --------------- RETURN RESULTS --------------- ! #
    
    
  } # R_grpnet_binomial.R



R_grpnet_ordinal_int <-
  function(nobs, nresp, eta, y, w, yfreq, theta, maxdif){
    
    # nobs = n
    # nresp = J
    # eta = x %*% beta
    # y = integer vector (elements in 1,...,nresp)
    # w = non-negative weights
    # yfreq = integer vector of length nresp
    # theta = shift vector of length J-1
    
    
    #########***#########  GRADIENT AND HESSIAN   #########***#########
    
    # update bigtheta 
    bigtheta <- c(100, theta, -100)
    
    # initialize grad, avec, and bvec
    grad <- rep(0.0, nresp-1L)
    avec <- rep(0.0, nresp-1L)
    bvec <- rep(0.0, nresp-2L)
    
    # define shift + linear predictor (for y and y+1)
    fit0 <- eta + bigtheta[y]
    fit1 <- eta + bigtheta[y+1L]
    
    # evaluate gradient / hessian components
    Pi0 <- 1 / (1 + exp(-fit0))
    Pi1 <- 1 / (1 + exp(-fit1))
    u0 <- Pi0 * (1 - Pi0)
    u1 <- Pi1 * (1 - Pi1)
    mu <- Pi0 - Pi1
    
    # update bvec (off-diagonals)
    dtop <- exp(bigtheta[1:nresp] + bigtheta[2:(nresp+1)])
    dbot <- ( exp(bigtheta[1:nresp]) - exp(bigtheta[2:(nresp+1)]) )^2
    dtheta <- dtop / dbot
    bvec <- (yfreq[2:(nresp-1)] * dtheta[2:(nresp-1)]) / nobs
    
    # update other shift gradient and hessian components
    for(i in 1:(nresp-1)){
      j <- i + 1L
      yj <- (y == j)
      yjm1 <- (y == (j-1))
      grad[i] <- mean( w * (u0 * yj - u1 * yjm1) / mu )
      vj <- 1 / (1 + exp(-(eta + bigtheta[j])))
      vj <- (1 - vj)^2 - (1 - vj)
      va <- vj - dtheta[j-1]
      vb <- vj - dtheta[j]
      avec[i] <- - mean(w * ( yjm1 * va + yj * vb ) )
    }
    
    
    #########***#########  INVSYMTRIDIAG   #########***#########
    
    # initializations
    n <- nresp - 1L
    u <- rep(0.0, n)
    v <- rep(0.0, n)
    d <- rep(0.0, n)
    delta <- rep(0.0, n)
    
    # update d
    d[n] <- avec[n]
    for(i in (n-1):1){
      d[i] <- avec[i] - bvec[i]^2 / d[i+1]
    }
    
    # update v
    v[1] <- 1 / d[1]
    for(i in 2:n){
      v[i] <- prod(bvec[1:(i-1)] / d[1:(i-1)]) / d[i]
    }
    
    # update delta
    delta[1] <- avec[1]
    for(i in 1:(n-1)){
      delta[i+1] <- avec[i+1] - bvec[i]^2 / delta[i]
    }
    
    # update u
    u[n] <- 1 / (delta[n] * v[n])
    for(i in 1:(n-1)){
      u[i] <- prod(bvec[i:(n-1)] / delta[i:(n-1)]) / (delta[n] * v[n])
    }
    
    
    #########***#########  NEWTON-RAPHSON UPDATE   #########***#########
    
    # calculate solve(hess) %*% grad
    ivec <- rep(0.0, n)
    for(i in 1:n){
      ivec[i] <- sum( (u[pmin(i, 1:n)] * v[pmax(i, 1:n)]) * grad )
    }
    
    # update maxdif and theta
    maxdif <- max( abs(ivec) / (1.0 + abs(theta)) )
    theta <- theta + ivec
    
    # return results
    res <- list(nobs = nobs,
                nresp = nresp,
                eta = eta, 
                y = y, 
                w = w,
                yfreq = yfreq,
                theta = theta,
                maxdif = maxdif)
    return(res)
    
  } # end R_grpnet_ordinal_int
