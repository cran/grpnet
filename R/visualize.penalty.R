visualize.penalty <- 
  function(x = seq(-5, 5, length.out = 1001), 
           penalty = c("LASSO", "MCP", "SCAD"), 
           alpha = 1, 
           lambda = 1, 
           gamma = 4, 
           derivative = FALSE,
           plot = TRUE,
           subtitle = TRUE,
           legend = TRUE,
           location = ifelse(derivative, "bottom", "top"),
           ...){
    # plot grpnet penalties
    # Nathaniel E. Helwig (helwig@umn.edu)
    # 2024-06-04
           
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    # check x
    xorig <- sort(as.numeric(x))
    x <- abs(xorig)
    nobs <- length(x)
    
    # check penalty
    penalties <- c("LASSO", "MCP", "SCAD")
    penalty <- sort(unique(as.character(penalty)))
    pid <- pmatch(penalty, penalties)
    if(any(is.na(pid))) stop(paste("Input 'penalty' is invalid. Must be some subset of\n", 
                                   paste(penalties, collapse = ", ")))
    penalty <- penalties[pid]
    npen <- length(penalty)
    
    # check alpha
    alpha <- as.numeric(alpha[1])
    if(alpha < 0 | alpha > 1) stop("Input 'alpha' must satisfy:  0 <= alpha <= 1")
    
    # check lambda
    lambda <- as.numeric(lambda[1])
    if(lambda < 0) stop("Input 'lambda' must satisfy:  lambda >= 0")
    
    # check gamma
    gamma <- as.numeric(gamma[1])
    if(any(penalty == "SCAD")){
      if(gamma <= 2) stop("Input 'gamma' must satisfy: gamma > 2 (SCAD)")
    } else if(any(penalty == "MCP")){
      if(gamma <= 1) stop("Input 'gamma' must satisfy: gamma > 1 (MCP)")
    } 
    
    # check derivative
    derivative <- as.logical(derivative[1])
    if(!any(derivative == c(TRUE, FALSE))) stop("Input 'derivative' must be TRUE/FALSE")
    
    
    
    #########***#########   EVALUATE PENALTY   #########***#########
    
    # initialize matrix to hold results
    res <- matrix(0.0, nrow = nobs, ncol = npen)
    colnames(res) <- penalty
    
    # loop through penalties
    for(j in 1:npen){
      
      if(derivative){
        if(penalty[j] == "LASSO"){
          res[,j] <- lambda
        } else if (penalty[j] == "MCP"){
          id <- (x <= gamma * lambda)
          res[id,j] <- lambda - x[id] / gamma
          res[!id,j] <- 0
        } else if (penalty[j] == "SCAD"){
          id <- which(x <= lambda)
          if(length(id) > 0) res[id,j] <- lambda
          id <- which((x > lambda) & (x <= gamma * lambda))
          if(length(id) > 0) res[id,j] <- (gamma * lambda - x[id]) / (gamma - 1)
          id <- which(x >= gamma * lambda)
          if(length(id) > 0) res[id,j] <- 0
        }
      } else {
        if(penalty[j] == "LASSO"){
          res[,j] <- lambda * x
        } else if (penalty[j] == "MCP"){
          id <- (x <= gamma * lambda)
          res[id,j] <- lambda * x[id] - x[id]^2 / (2 * gamma)
          res[!id,j] <- gamma * lambda^2 / 2
        } else if (penalty[j] == "SCAD"){
          id <- which(x <= lambda)
          if(length(id) > 0) res[id,j] <- lambda * x[id]
          id <- which((x > lambda) & (x <= gamma * lambda))
          if(length(id) > 0) res[id,j] <- (gamma * lambda * x[id] - 0.5 * (x[id]^2 + lambda^2)) / (gamma - 1)
          id <- which(x >= gamma * lambda)
          if(length(id) > 0) res[id,j] <- lambda^2 * (gamma + 1) / 2
        }
      } # end if(derivative)
      
      # add ridge component?
      if(alpha < 1){
        if(derivative){
          res <- alpha * res + (1-alpha) * abs(x)
        } else {
          res <- alpha * res + (1-alpha) * x^2 / 2
        }
      }
      
    } # end for(j in 1:npen)
    
    # return results?
    if(!plot){
      return(as.data.frame(cbind(x = xorig, res))) 
    }
    
    
    
    #########***#########   COLLECT ARGS   #########***#########
      
    # collect ...
    args <- list(...)
    
    # check args$xlim
    if(is.null(args$xlim)) args$xlim <- extendrange(xorig)
    
    # check args$ylim
    if(is.null(args$ylim)) args$ylim <- extendrange(res)
    
    # check args$lty
    if(is.null(args$lty)) args$lty <- c(1, 2, 4)
    
    # check args$lwd
    if(is.null(args$lwd)) args$lwd <- c(2, 2, 2)
    
    # check args$col
    if(is.null(args$col)) args$col <- c("black", "blue", "red")
    
    # check args$xlab
    if(is.null(args$xlab)) args$xlab <- expression(theta)
    
    # check args$ylab
    if(is.null(args$ylab)) args$ylab <- ifelse(derivative,
                                               expression(italic(P)*"'" * (theta)),
                                               expression(italic(P)(theta)))
    # check args$main
    if(is.null(args$main)) args$main <- ifelse(derivative, "penalty derivatives", "penalty functions")
    
    
    
    #########***#########   PLOT PENALTY   #########***#########
    
    # plot first penalty
    plot(xorig, res[,1], type = "l", lty = args$lty[1], col = args$col[1], 
         lwd = args$lwd[1], xlim = args$xlim, ylim = args$ylim, 
         xlab = args$xlab, ylab = args$ylab, main = args$main)
    
    if(npen > 1L){
      for(j in 2:npen){
        lines(xorig, res[,j], lty = args$lty[j], col = args$col[j], lwd = args$lwd[j])
      }
    }
    
    # subtitle?
    if(subtitle){
      subtit <- bquote(alpha == .(alpha) ~ "    " ~ gamma == .(gamma) ~ "    " ~ lambda == .(lambda))
      mtext(as.expression(subtit), side = 3)
    }
    
    # legend?
    if(legend){
      legend(location, legend = penalty, 
             col = args$col, lty = args$lty, lwd = args$lwd, bty = "n")
    }
    
    
  } # end visualize.penalty