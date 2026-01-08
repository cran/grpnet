family.grpnet <-
  function(object, theta = 1){
    # prepare family for grpnet
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2025-08-15
    
    # check class of object
    if(inherits(object, "cv.grpnet")){
      return(object$grpnet.fit$family)
    } else if(inherits(object, "grpnet")){
      return(object$family)
    }
    
    # check input family object
    object <- as.character(object[1])
    if(object == "mgaussian" | object == "mvn") object <- "multigaussian"
    if(object == "hsvm") object <- "svm1"
    if(object == "sqsvm") object <- "svm2"
    if(object == "polr") object <- "ordinal"
    families <- c("gaussian", "multigaussian", "svm1", "svm2", "logit", "binomial", "multinomial", "ordinal", "poisson", "negative.binomial", "Gamma", "inverse.gaussian")
    object <- pmatch(object, families)
    if(is.na(object)) stop("'object' not recognized")
    object <- families[object]
    
    # create needed family info
    if(object %in% c("gaussian", "multigaussian")){
      
      famname <- object
      object <- gaussian()
      object$family <- famname
      
    } else if(object == "svm1"){
      
      theta <- as.numeric(theta[1])
      if(theta <= 0) stop("Input 'theta' must be positive")
      .Theta <- theta
      env <- new.env(parent = .GlobalEnv)
      assign(".Theta", theta, envir = env)
      dr <- function(y, mu, wt){
        muy <- mu * y
        vec <- rep(0.0, length(muy))
        idx <- which(muy > 1.0 - .Theta)
        vec[idx] <- (1.0 - muy[idx])^2 / (2.0 * .Theta)
        vec[!idx] <- 1.0 - muy[!idx] - .Theta / 2.0
        vec[which(muy > 1.0)] <- 0.0
        2 * wt * vec
      }
      object <- list(family = "svm1",
                     linkinv = function(eta) eta,
                     dev.resids = dr)
      
    } else if(object == "svm2"){
      
      dr <- function(y, mu, wt){
        2 * wt * pmax(0, 1 - mu * y)^2
      }
      object <- list(family = "svm2",
                     linkinv = function(eta) eta,
                     dev.resids = dr)
      
    } else if(object == "logit"){
      
      dr <- function(y, mu, wt){
        2 * wt * log(1 + exp(-mu*y)) 
      }
      object <- list(family = "logit",
                     linkinv = function(eta) {1 / (1 + exp(-eta))},
                     dev.resids = dr)
      
    } else if(object == "binomial"){
      
      dr <- function(y, mu, wt){
        ylogy <- function(y, mu){
          res <- y * log(y / mu)
          res[is.nan(res)] <- 0
          res
        }
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        2 * wt * (ylogy(y, mu) + ylogy(1 - y, 1 - mu))
      }
      object <- list(family = "binomial",
                     linkinv = function(eta) {1 / (1 + exp(-eta))},
                     dev.resids = dr)
      
    } else if (object == "multinomial"){
      
      il <- function(eta){
        expeta <- exp(eta - apply(eta, 1, max))
        mu <- expeta / rowSums(expeta)
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        mu
      } # end il
      dr <- function(y, mu, wt) {
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        -2 * wt * rowSums(y * log(mu))
      }
      object <- list(family = "multinomial",
                     linkinv = il,
                     dev.resids = dr)
      
    } else if(object == "ordinal"){
      
      dr <- function(y, mu, wt){
        mu[mu < 0.000001] <- 0.000001
        mu[mu > 0.999999] <- 0.999999
        -2 * wt * log(mu)
      }
      object <- list(family = "ordinal",
                     linkinv = function(eta) {1 / (1 + exp(-eta))},
                     dev.resids = dr)
      
    } else if(object == "poisson"){
      
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        r <- mu * wt
        p <- which(y > 0)
        if(is.matrix(mu) && ncol(mu) > 1L){
          r[p,] <- (wt * (y * log(y/mu) - (y - mu)))[p,]
        } else {
          r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
        }
        2 * r
      }
      object <- list(family = "poisson",
                     linkinv = il,
                     dev.resids = dr)
      
    } else if(object == "negative.binomial"){
      
      theta <- as.numeric(theta[1])
      if(theta <= 0) stop("Input 'theta' must be positive")
      .Theta <- theta
      env <- new.env(parent = .GlobalEnv)
      assign(".Theta", theta, envir = env)
      
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        2 * wt * ( y * log(pmax(1, y) / mu) - (y + .Theta) * log((y + .Theta) / (mu + .Theta)) )
      }
      object <- list(family = "negative.binomial",
                     linkinv = il,
                     dev.resids = dr)
      
    } else if(object == "Gamma"){
      
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        -2 * wt * (log(y / mu) - (y - mu) / mu)
      }
      object <- list(family = "Gamma",
                     linkinv = il,
                     dev.resids = dr)
      
    } else if(object == "inverse.gaussian"){
      
      il <- function(eta) pmax(exp(eta), .Machine$double.eps)
      dr <- function(y, mu, wt){
        wt * ( (y - mu)^2 / (y * mu^2) )
      }
      object <- list(family = "inverse.gaussian",
                     linkinv = il,
                     dev.resids = dr)
      
    } # end if(object == "gaussian")
    
    return(object)
    
  } # family.grpnet.R