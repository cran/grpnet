# Some code (for interpolating predictions) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

predict.grpnet <-
  function(object, 
           newx,
           newdata,
           s = NULL,
           type = c("link", "response", "class"),
           ...){
    # predict from a fit grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-07-05
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### check object
    if(!inherits(object, "grpnet")) stop("Input 'object' must be of class 'grpnet'")
  
    ### check newx / newdata
    if(is.null(object$formula)){
      newx <- as.matrix(newx)
      ncoefs <- if(is.list(object$beta)) nrow(object$beta[[1]]) else nrow(object$beta)
      if(ncol(newx) != ncoefs) stop("Input 'newx' must satisfy:  ncol(newx) == nrow(object$beta)")
    } else {
      newx <- model.matrix(object = object$formula[-2], data = newdata)
      if(colnames(newx)[1] == "(Intercept)") newx <- newx[,-1]
      ncoefs <- if(is.list(object$beta)) nrow(as.matrix(object$beta[[1]])) else nrow(as.matrix(object$beta))
      if(ncol(newx) != ncoefs) stop("Input 'newdata' produced a design matrix of the wrong dimension\n(likely due to a factor level mismatch between 'data' and 'newdata')")
    }
    nobs <- nrow(newx)
    newxnames <- rownames(newx)
    if(is.null(newxnames)) newxnames <- 1:nobs
    
    ### check s
    if(is.null(s)){
      newlambdas <- FALSE
      s <- object$lambda
    } else {
      newlambdas <- TRUE
      s <- as.numeric(s)
      if(any(s < 0)) stop("Inpust 's' must be a vector of non-negative numerics")
    } # end if(is.null(s))
    
    ### number of fit lambdas
    nlam <- length(object$lambda)
    
    ### number of new lambdas
    ns <- length(s)
    
    ### check type
    family <- object$family$family
    type <- pmatch(as.character(type[1]), c("link", "response", "class"))
    if(is.na(type)) stop("Invalid 'type' input")
    type <- c("link", "response", "class")[type]
    if(type == "class" && !(family %in% c("binomial", "multinomial")))
      stop("Input 'type' can only be set to 'class' for binomial and multinomial families")
    
    
    
    ######***######   MULTINOMIAL   ######***######
    
    if(family == "multinomial"){ 
      
      # note: needs special code to return predictions for each response class
      
      ### get number of response classes
      nresp <- length(object$beta)
      
      ### initialize array for predictions
      fit <- array(data = NA, dim = c(nobs, nresp, ns))
      dimnames(fit) <- list(newxnames, object$ylev, paste0("s", 1:ns))
      
      ### predictions at object$lambda are easy...
      if(!newlambdas){
        for(k in 1:nresp) {
          if(nlam == 1L){
            fit[,k,] <- cbind(1, newx) %*% c(object$a0[k], object$beta[[k]])
          } else {
            fit[,k,] <- cbind(1, newx) %*% rbind(object$a0[k,], object$beta[[k]])
          }
        }
        if(type %in% c("response", "class")) {
          for(i in 1:ns) fit[,,i] <- object$family$linkinv(fit[,,i])
          if(type == "class") fit <- matrix(object$ylev[apply(fit, c(1,3), which.max)], nrow = nobs, ncol = ns)
        }
        return(drop(fit))
      }
      
      ### only 1 fit lambda?
      if(nlam == 1L){
        for(k in 1:nresp){
          fit[,k,] <- cbind(1, newx) %*% c(object$a0[k], object$beta[[k]])
        }
        if(type == "response") {
          fit <- object$family$linkinv(fit)
        } else if(type == "class"){
          # note: family == "multinomial" is implied
          fit <- object$ylev[apply(object$family$linkinv(fit), 1, which.max)]
        }
        return(fit)
      }
      
      ### min and max lambda from fit model
      lambda.max <- max(object$lambda)
      lambda.min <- min(object$lambda)
      
      ### transform lambda (and s) to the [0,1] interval
      sfrac <- (lambda.max - s) / (lambda.max - lambda.min)
      lambda <- (lambda.max - object$lambda) / (lambda.max - lambda.min)
      
      ### correct for out-of-range s values
      sfrac[sfrac > max(lambda)] <- max(lambda)
      sfrac[sfrac < min(lambda)] <- min(lambda)
      
      ### linear interpolation at sfrac
      interp <- approx(x = lambda, y = seq(lambda), xout = sfrac)$y
      
      ### L and R indices
      left <- floor(interp)
      right <- ceiling(interp)
      
      ### update sfrac w.r.t. L and R indices
      sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
      sfrac[left == right] <- 1
      sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
      
      ### interpolate...
      for(k in 1:nresp){
        coefs <- rbind(object$a0[k,], object$beta[[k]])
        coefs <- coefs[, left, drop=FALSE] %*% diag(sfrac, nrow = ns, ncol = ns) + coefs[, right, drop=FALSE] %*% diag(1 - sfrac, nrow = ns, ncol = ns)
        colnames(coefs) <- paste0("s", 1:ns)
        fit[,k,] <- cbind(1, newx) %*% coefs
      }
      
      ### return fitted values
      if(type %in% c("response", "class")) {
        for(i in 1:ns) fit[,,i] <- object$family$linkinv(fit[,,i])
        if(type == "class") fit <- matrix(object$ylev[apply(fit, c(1,3), which.max)], nrow = nobs, ncol = ns)
      }
      return(drop(fit))
      
    } # end if(object$family$family == "multinomial")
    
    
    
    ######***######   OTHER FAMILIES   ######***######
    
    ### predictions at object$lambda are easy...
    if(!newlambdas){
      if(nlam == 1L){
        fit <- cbind(1, newx) %*% c(object$a0, object$beta)
      } else {
        fit <- cbind(1, newx) %*% rbind(object$a0, object$beta)
      }
      if(type == "response") {
        fit <- object$family$linkinv(fit)
      } else if(type == "class"){ 
        # note: family == "binomial" is implied
        fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
      }
      return(drop(fit))
    }
    
    ### only 1 fit lambda?
    if(nlam == 1L){
      fit <- cbind(1, newx) %*% c(object$a0, object$beta)
      if(type == "response") {
        fit <- object$family$linkinv(fit)
      } else if(type == "class"){
        # note: family == "binomial" is implied
        fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
      }
      return(drop(matrix(fit, nrow = nrow(newx), ncol = ns)))
    }
    
    ### min and max lambda from fit model
    lambda.max <- max(object$lambda)
    lambda.min <- min(object$lambda)
    
    ### transform lambda (and s) to the [0,1] interval
    sfrac <- (lambda.max - s) / (lambda.max - lambda.min)
    lambda <- (lambda.max - object$lambda) / (lambda.max - lambda.min)
    
    ### correct for out-of-range s values
    sfrac[sfrac > max(lambda)] <- max(lambda)
    sfrac[sfrac < min(lambda)] <- min(lambda)
    
    ### linear interpolation at sfrac
    interp <- approx(x = lambda, y = seq(lambda), xout = sfrac)$y
    
    ### L and R indices
    left <- floor(interp)
    right <- ceiling(interp)
    
    ### update sfrac w.r.t. L and R indices
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
    
    ### interpolate...
    coefs <- rbind(object$a0, object$beta)
    coefs <- coefs[, left, drop=FALSE] %*% diag(sfrac, nrow = ns, ncol = ns) + coefs[, right, drop=FALSE] %*% diag(1 - sfrac, nrow = ns, ncol = ns)
    colnames(coefs) <- paste0("s", 1:ns)
    
    ### return fitted values
    fit <- cbind(1, newx) %*% coefs
    if(type == "response") {
      fit <- object$family$linkinv(fit)
    } else if(type == "class"){
      # note: family == "binomial" is implied
      fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
    }
    return(drop(fit))
    
  } # end predict.grpnet