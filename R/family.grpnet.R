family.grpnet <-
  function(object, theta = 1){
    # prepare family for grpnet
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2024-06-28
    
    object <- as.character(object[1])
    if(object == "gaussian"){
      
      object <- gaussian()
      
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