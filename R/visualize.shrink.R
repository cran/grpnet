visualize.shrink <-
  function(x = seq(-5, 5, length.out = 1001), 
           penalty = c("LASSO", "MCP", "SCAD"), 
           alpha = 1, 
           lambda = 1, 
           gamma = 4, 
           fitted = FALSE,
           plot = TRUE,
           subtitle = TRUE,
           legend = TRUE,
           location = "top",
           ...){
    # plot grpnet shrinkage and selection operators
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
    
    
    
    #########***#########   EVALUATE SHRINKAGE   #########***#########
    
    # initialize l1 and l2 penalty weights
    lambda1 <- lambda * alpha
    lambda2 <- lambda * (1-alpha)
    
    # initialize matrix to hold results
    res <- matrix(0.0, nrow = nobs, ncol = npen)
    colnames(res) <- penalty
    
    # loop through penalties
    for(j in 1:npen){
      
      if(penalty[j] == "LASSO"){
        res[,j] <- pmax(0, 1 - lambda1 / x) / (1 + lambda2)
      } else if (penalty[j] == "MCP"){
        id <- (x <= gamma * lambda1 * (1 + lambda2))
        res[id,j] <- pmax(0, 1 - lambda1 / x[id]) / (1 + lambda2 - 1/gamma)
        res[!id,j] <- 1 / (1 + lambda2)
      } else if (penalty[j] == "SCAD"){
        res[,j] <- 1 / (1 + lambda2)
        id <- which(x <= (lambda1 + lambda1 * (1 + lambda2)))
        res[id,j] <- pmax(0, 1 - lambda1 / x[id]) / (1 + lambda2)
        id <- which( (x > (lambda1 + lambda1 * (1 + lambda2))) & (x <= gamma * lambda1 * (1 + lambda2)) )
        res[id,j] <- pmax(0, 1 - (lambda1 / x[id]) * (gamma/(gamma-1)) ) / (1 + lambda2 - 1 / (gamma-1))
      }
      
    }  # end for(j in 1:npen)
    
    # fitted values instead of shrinkage factors?
    if(fitted) res <- xorig * res
    
    # return results?
    if(!plot){
      return(as.data.frame(cbind(x = xorig, res))) 
    }
    
    
    
    #########***#########   COLLECT ARGS   #########***#########
    
    # collect ...
    args <- list(...)
    
    # check args$ylim
    if(is.null(args$ylim)) args$ylim <- extendrange(res)
    
    # check args$xlim
    if(is.null(args$xlim) && fitted) args$xlim <- args$ylim
    
    # check args$lty
    if(is.null(args$lty)) args$lty <- c(1, 2, 4)
    
    # check args$lwd
    if(is.null(args$lwd)) args$lwd <- c(2, 2, 2)
    
    # check args$col
    if(is.null(args$col)) args$col <- c("black", "blue", "red")
    
    # check args$xlab
    if(is.null(args$xlab)) args$xlab <- expression(theta)
    
    # check args$ylab
    if(is.null(args$ylab)) args$ylab <- ifelse(fitted,
                                               expression(italic(S)[lambda[1] * ", " * lambda[2]](theta) %*% theta),
                                               expression(italic(S)[lambda[1] * ", " * lambda[2]](theta)))
      
    # check args$main
    if(is.null(args$main)) args$main <- ifelse(fitted,
                                               "shrunken estimators",
                                               "shrinkage operators")
    
    
    
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
    
    # fitted?
    if(fitted) abline(a = 0, b = 1, lty = 3)
    
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
    
    
  } # end visualize.shrink