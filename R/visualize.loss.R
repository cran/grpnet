visualize.loss <-
  function(x = seq(-3, 3, length.out = 1001),
           family = c("gaussian", "multigaussian",
                      "svm1", "svm2", "logit",
                      "binomial", "multinomial", 
                      "poisson", "negative.binomial", 
                      "Gamma", "inverse.gaussian"),
           theta = 1,
           type = c("link", "response"),
           y = NULL,
           plot = TRUE,
           add = FALSE,
           ...){
    # plot grpnet loss functions
    # Nathaniel E. Helwig (helwig@umn.edu)
    # 2025-06-03
    
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    # check family
    family <- as.character(family[1])
    fam <- family.grpnet(family, theta = theta)
    
    # check theta
    theta <- as.numeric(theta[1])
    if(theta <= 0.0) stop("Input 'theta' must be positive")
    
    # check type
    type <- as.character(type[1])
    types <- c("link", "response")
    type <- pmatch(type, types)
    if(is.na(type)) stop("Input 'type' must be 'link' or 'response'")
    type <- types[type]
    
    # check x and y
    x <- as.numeric(x)
    if(!is.null(y)) y <- as.numeric(y[1])
    
    # check add
    add <- as.logical(add[1])
    if(!any(add == c(TRUE, FALSE))) stop ("Input 'add' must be TRUE or FALSE")
    
    # convert x to mu
    mu <- fam$linkinv(x)
    
    
    
    #########***#########   EVALUATE LOSS   #########***#########
    
    # evaluate loss
    if(family %in% c("gaussian", "multigaussian")){
      y <- ifelse(is.null(y), 0, y[1])
      loss <- (y - mu)^2
    } else if(family == "svm1"){
      y <- ifelse(is.null(y), 1, y[1])
      muy <- mu * y
      loss <- rep(0, length(mu))
      id <- (muy > 1 - theta)
      loss[id] <- pmax(1 - muy[id], 0)^2 / (2 * theta)
      loss[!id] <- 1 - muy[!id] - theta/2
    } else if(family == "svm2"){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- pmax(0, 1 - mu * y)^2
    } else if(family == "logit"){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- log(1 + exp(-x * y))
    } else if(family %in% c("binomial", "multinomial")){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- - y * log(mu) - (1 - y) * log(1 - mu)
    } else if(family == "poisson"){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- mu - y * log(mu)
    } else if(family == "negative.binomial"){
      y <- ifelse(is.null(y), 1, y[1])
      const <- lgamma(theta) - lgamma(theta + y) - theta * log(theta)
      loss <- (theta + y) * log(theta + mu) - y * log(mu) + const
    } else if(family == "Gamma"){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- log(mu) + y / mu
    } else if(family == "inverse.gaussian"){
      y <- ifelse(is.null(y), 1, y[1])
      loss <- (y - mu)^2 / (mu^2 * y)
    }
    
    
    
    #########***#########   RETURN LOSS?   #########***#########
    
    if(!plot){
      df <- data.frame(eta = x, mu = mu, loss = loss)
      return(df)
    } 
    
    
    
    #########***#########   PLOT LOSS   #########***#########
    
    # collect ...
    args <- list(...)
    
    # add x and y
    if(type == "link"){
      args$x <- x
    } else {
      args$x <- mu
    }
    args$y <- loss
    
    # add xlab and ylab
    if(is.null(args$xlab)) {
      if(type == "link"){
        args$xlab <- expression(italic(eta))
      } else {
        args$xlab <- expression(italic(mu))
      }
    }
    if(is.null(args$ylab)) {
      if(type == "link"){
        args$ylab <- substitute(expression(italic(L) * "( " * eta * " | " * italic(y) * " = " * yval * " )"),
                                list(yval = y))
      } else {
        args$ylab <- substitute(expression(italic(L) * "( " * mu * " | " * italic(y) * " = " * yval * " )"),
                                list(yval = y))
      }
    }
    if(is.null(args$main)) args$main <- family
    
    # add type
    if(!add) args$type <- "l"
    
    # check args$xlim
    if(is.null(args$xlim)) {
      if(type == "link"){
        args$xlim <- extendrange(x)
      } else {
        args$xlim <- extendrange(mu)
      }
    }
    
    # check args$ylim
    if(is.null(args$ylim)) args$ylim <- extendrange(loss)
    
    # check args$lty
    if(is.null(args$lty)) args$lty <- 1L
    
    # check args$lwd
    if(is.null(args$lwd)) args$lwd <- 2L
    
    # check args$col
    if(is.null(args$col)) args$col <- "darkgray"
    
    # add lines or draw plot
    if(add){
      do.call(lines, args)
    } else {
      rm(plot)  # remove logical "plot" argument
      do.call(plot, args)
    }
    
  } # end visualize.loss