R_grpnet_multinom <-
  function(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw, 
           alpha, nlam, lambda, lmr, penid, gamma, eps, maxit,
           standardize, intercept, ibeta, betas, iters,
           nzgrps, nzcoef, edfs, devs, nulldev){
    # grpnet_multinom.f90 translation to R
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2024-06-28
    
    
    # ! --------------- LOCAL DEFINITIONS --------------- ! #
    ia <- ib <- rep(NA, ngrps)
    xmean <- rep(NA, nvars)
    xev <- rep(NA, ngrps)
    difbeta <- matrix(NA, nvars, nresp)
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
    xev <- xev / 4
    # ! --------------- GET MAX EIGENVALUE --------------- ! #
    
    
    # ! --------------- MISCELLANEOUS INITIALIZATIONS --------------- ! #
    maxlam <- max(lambda)
    makelambda <- 0L
    iter <- 0L
    active <- strong <- rep(0L, ngrps)
    nzgrps <- nzcoef <- rep(0L, nlam)
    ibeta <- matrix(0.0, nresp, nlam)
    beta <- zvec <- grad <- matrix(0.0, nvars, nresp)
    gradnorm <- rep(0.0, ngrps)
    twolam <- 0.0
    devs <- rep(0.0, nlam)
    eta <- off
    expeta <- exp( eta - apply(eta, 1, max) )
    mu <- expeta / rowSums(expeta)
    r <- w * (y - mu)
    # ! --------------- MISCELLANEOUS INITIALIZATIONS --------------- ! #
    
    
    # ! --------------- GENERATE LAMBDA --------------- ! #
    if(maxlam <= macheps){
      
      makelambda <- 1L
      i <- 1L
      
      # ! find unpenalized groups ! # 
      for(k in 1:ngrps){
        if(pw[k] <= macheps){
          active[k] <- 1L
          nzgrps[i] <- nzgrps[i] + 1L
          nzcoef[i] <- nzcoef[i] + gsize[k] * nresp
        }
      }
      # ! find unpenalized groups ! #
      
      # ! iterate until active coefficients converge ! #
      if(nzgrps[i] > 0L){
        
        while(iter < maxit){
          
          # ! update iter and reset counters
          ctol <- 0.0
          iter <- iter + 1L
          
          # ! update active groups
          for(k in 1:ngrps){
            if(active[k] == 0L) next
            grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
            zvec[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + grad[ia[k]:ib[k],] / xev[k]
            zvec[ia[k]:ib[k],] <- zvec[ia[k]:ib[k],] - rowMeans(zvec[ia[k]:ib[k],,drop=FALSE])
            difbeta[ia[k]:ib[k],] <- zvec[ia[k]:ib[k],] - beta[ia[k]:ib[k],]
            maxdif <- max( abs(difbeta[ia[k]:ib[k],]) / (1.0 + abs(difbeta[ia[k]:ib[k],])) )
            beta[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + difbeta[ia[k]:ib[k],]
            eta <- eta + (x[,ia[k]:ib[k]] %*% difbeta[ia[k]:ib[k],]) / w
            expeta <- exp( eta - apply(eta, 1, max) )
            mu <- expeta / rowSums(expeta)
            r <- w * (y - mu)
            ctol <- max(maxdif, ctol)
          } # ! k=1,ngrps
          
          # ! update intercept
          if(intercept == 1L){
            zint <- ibeta[,i] + ( colSums(r * w) / nobs ) * 4
            zint <- zint - sum(zint) / nresp
            difibeta <- zint - ibeta[,i]
            maxdif <- abs(difibeta) / (1.0 + abs(ibeta[,i]))
            ibeta[,i] <- ibeta[,i] + difibeta
            eta <- eta + matrix(difibeta, nobs, nresp, byrow = TRUE)
            expeta <- exp( eta - apply(eta, 1, max) )
            mu <- expeta / rowSums(expeta)
            r <- w * (y - mu)
            ctol <- max(maxdif, ctol)
          } # if(intercept == 1L)
          
          # ! convergence check
          if(ctol < eps) break
          
        } # ! WHILE(iter < maxit)
        
      } else {
        
        # ! intercept only
        iter <- 1L
        if(intercept == 1L){
          ibeta[,i] <- log( colSums(y * w^2) / nobs )
          ibeta[,i] <- ibeta[,i] - sum(ibeta[,i]) / nresp
          eta <- eta + matrix(ibeta[,i], nobs, nresp, byrow = TRUE)
          expeta <- exp( eta - apply(eta, 1, max) )
          mu <- expeta / rowSums(expeta)
          r <- w * (y - mu)
        }
        
      } # ! (nzgrps[i] > 0)
      # ! iterate until active coefficients converge ! #
      
      # ! create lambda sequence ! #
      for(k in 1:ngrps){
        if(pw[k] > macheps){
          grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
          gradnorm[k] <- sqrt( sum(grad[ia[k]:ib[k],]^2) ) / pw[k]
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
      devs[i] <- R_grpnet_multinom_dev(nobs, nresp, y, mu, w^2)
      # ! calculate deviance ! #
      
      # ! save results ! # 
      betas[,,i] <- beta
      iters[i] <- iter
      nzgrps[i] <- nzgrps[i] + intercept
      nzcoef[i] <- nzcoef[i] + intercept * nresp
      edfs[i] <- intercept * (nresp - 1)
      # ! save results ! # 
      
    }
    # ! --------------- GENERATE LAMBDA --------------- ! #
    
    
    # ! --------------- ITERATIVE WORK --------------- ! #
    for(i in 1:nlam){
      
      # ! initializations ! # 
      if(i == 1L && makelambda == 1L) next
      if(i > 1L){
        ibeta[,i] <- ibeta[,i-1]
        beta <- betas[,,i-1]
        twolam <- alpha * (2.0 * lambda[i] - lambda[i-1])
      } else {
        grad <- crossprod(x, r) / nobs
      }
      # ! initializations ! # 
      
      # ! strong rule initialization ! #
      for(k in 1:ngrps){
        gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k],]^2))
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
            
            # ! unweighted or weighted update?
            if(weighted == 0L){
              
              # ! update active groups
              for(k in 1:ngrps){
                if(active[k] == 0L) next
                penone <- alpha * lambda[i] * pw[k] / xev[k]
                pentwo <- (1.0 - alpha) * lambda[i] * pw[k] / xev[k]
                grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
                zvec[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + grad[ia[k]:ib[k],] / xev[k]
                zvec[ia[k]:ib[k],] <- zvec[ia[k]:ib[k],] - rowMeans(zvec[ia[k]:ib[k],,drop=FALSE])
                znorm <- sqrt(sum(zvec[ia[k]:ib[k],]^2))
                bnorm <- sqrt(sum(beta[ia[k]:ib[k],]^2))
                shrink <- R_grpnet_penalty(znorm, penid, penone, pentwo, gamma)
                if(shrink == 0.0 && bnorm == 0.0) next
                difbeta[ia[k]:ib[k],] <- shrink * zvec[ia[k]:ib[k],] - beta[ia[k]:ib[k],]
                maxdif <- max( abs(difbeta[ia[k]:ib[k],]) / (1.0 + abs(beta[ia[k]:ib[k],])) )
                beta[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + difbeta[ia[k]:ib[k],]
                eta <- eta + x[,ia[k]:ib[k]] %*% difbeta[ia[k]:ib[k],]
                expeta <- exp( eta - apply(eta, 1, max) )
                mu <- expeta / rowSums(expeta)
                r <- y - mu
                ctol <- max(maxdif , ctol)
                if(shrink > 0.0){
                  nzgrps[i] <- nzgrps[i] + 1L
                  edfs[i] <- edfs[i] + gsize[k] * (nresp - 1) * shrink
                }
              } # ! k=1,ngrps
              
              # ! update intercept
              if(intercept == 1L){
                zint <- ibeta[,i] + ( colSums(r) / nobs ) * 4
                zint <- zint - sum(zint) / nresp
                difibeta <- zint - ibeta[,i]
                maxdif <- abs(difibeta) / (1 + abs(ibeta[,i]))
                ibeta[,i] <- ibeta[,i] + difibeta
                eta <- eta + matrix(difibeta, nobs, nresp, byrow = TRUE)
                expeta <- exp( eta - apply(eta, 1, max) )
                mu <- expeta / rowSums(expeta)
                r <- y - mu
                ctol <- max(maxdif, ctol)
              } # ! (intercept == 1)
              
            } else {
              
              # ! update active groups
              for(k in 1:ngrps){
                if(active[k] == 0L) next
                penone <- alpha * lambda[i] * pw[k] / xev[k]
                pentwo <- (1 - alpha) * lambda[i] * pw[k] / xev[k]
                grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
                zvec[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + grad[ia[k]:ib[k],] / xev[k]
                zvec[ia[k]:ib[k],] <- zvec[ia[k]:ib[k],] - rowMeans(zvec[ia[k]:ib[k],,drop=FALSE])
                znorm <- sqrt(sum(zvec[ia[k]:ib[k],]^2))
                bnorm <- sqrt(sum(beta[ia[k]:ib[k],]^2))
                shrink <- R_grpnet_penalty(znorm, penid, penone, pentwo, gamma)
                if(shrink == 0.0 && bnorm == 0.0) next
                difbeta[ia[k]:ib[k],] <- shrink * zvec[ia[k]:ib[k],] - beta[ia[k]:ib[k],]
                maxdif <- max( abs(difbeta[ia[k]:ib[k],]) / (1.0 + abs(beta[ia[k]:ib[k],])) )
                beta[ia[k]:ib[k],] <- beta[ia[k]:ib[k],] + difbeta[ia[k]:ib[k],]
                eta <- eta + (x[,ia[k]:ib[k]] %*% difbeta[ia[k]:ib[k],]) / w
                expeta <- exp( eta - apply(eta, 1, max) )
                mu <- expeta / rowSums(expeta)
                r <- w * (y - mu)
                ctol <- max(maxdif , ctol)
                if(shrink > 0.0){
                  nzgrps[i] <- nzgrps[i] + 1L
                  edfs[i] <- edfs[i] + gsize[k] * (nresp - 1) * shrink
                }
              } # ! k=1,ngrps
              
              # ! update intercept
              if(intercept == 1L){
                zint <- ibeta[,i] + ( colSums(r * w) / nobs ) * 4
                zint <- zint - sum(zint) / nresp
                difibeta <- zint - ibeta[,i]
                maxdif <- abs(difibeta) / (1 + abs(ibeta[,i]))
                ibeta[,i] <- ibeta[,i] + difibeta
                eta <- eta + matrix(difibeta, nobs, nresp, byrow = TRUE)
                expeta <- exp( eta - apply(eta, 1, max) )
                mu <- expeta / rowSums(expeta)
                r <- w * (y - mu)
                ctol <- max(maxdif, ctol)
              } # ! (intercept == 1)
              
            } # ! IF (weighted == 0)
            
            # ! convergence check
            if(ctol < eps) break
            
          } # ! WHILE(iter < maxit) - inner
          # ! iterate until active coefficients converge ! #
          
          # ! check inactive groups in strong set ! #
          violations <- 0L
          for(k in 1:ngrps){
            if(strong[k] == 0L | active[k] == 1L) next
            grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
            gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k],]^2))
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
          grad[ia[k]:ib[k],] <- crossprod(x[,ia[k]:ib[k]], r) / nobs
          gradnorm[k] <- sqrt(sum(grad[ia[k]:ib[k],]^2))
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
        for(j in 1:nresp){
          for(l in 1:gsize[k]){
            if(abs(beta[ia[k] + l - 1,j]) > macheps){
              nzcoef[i] <- nzcoef[i] + 1L
            }
          }
        }
      }
      # ! calculate nzcoef ! #
      
      # ! calculate deviance ! #
      devs[i] <- R_grpnet_multinom_dev(nobs, nresp, y, mu, w^2)
      # ! calculate deviance ! #
      
      # ! save results ! #
      betas[,,i] <- beta
      iters[i] <- iter
      nzgrps[i] <- nzgrps[i] + intercept
      nzcoef[i] <- nzcoef[i] + intercept * nresp
      edfs[i] <- edfs[i] + intercept * (nresp - 1)
      # ! save results ! #
      
    }
    # ! --------------- ITERATIVE WORK --------------- ! #
    
    
    # ! --------------- POST PROCESSING --------------- ! #
    if(standardize == 1L){
      for(k in 1:ngrps){
        betas[ia[k]:ib[k],,] <- betas[ia[k]:ib[k],,] / xsdev[k]
      }
    }
    if(intercept == 1L){
      for(j in 1:nresp){
        ibeta[j,] <- ibeta[j,] - xmean %*% betas[,j,]
      }
      mu <- matrix(colSums(y * w^2) / nobs, nobs, nresp, byrow = TRUE)
    } else {
      expeta <- exp( off - apply(off, 1, max) )
      mu <- expeta / rowSums(expeta)
    }
    nulldev <- R_grpnet_multinom_dev(nobs, nresp, y, mu, w^2)
    names(xsdev) <- names(pw)
    pw <- xsdev
    # ! --------------- POST PROCESSING --------------- ! #
    
    
    # ! --------------- RETURN RESULTS --------------- ! #
    res <- list(nobs = nobs, 
                nvars = nvars,
                nresp = nresp,
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
    
    
  } # R_grpnet_multinom.R


R_grpnet_multinom_dev <-
  function(nobs, nresp, y, mu, wt, dev){
    dev <- 0.0
    for(i in 1:nobs){
      for(j in 1:nresp){
        if(mu[i,j] < 1e-6){
          mu[i,j] <- 1e-6
        } else if(mu[i,j] > (1 - 1e-6)){
          mu[i,j] <- 1 - 1e-6
        }
      }
      dev <- dev - 2 * wt[i] * sum( y[i,] * log(mu[i,]) )
    }
    return(dev)
  } # R_grpnet_multinom_dev.R