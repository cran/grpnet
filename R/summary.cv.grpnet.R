summary.cv.grpnet <-
  function(object, ...){
    # summarize cv.grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2026-04-28
    
    if(!inherits(object, "cv.grpnet")) stop("Input 'object' should be of class 'cv.grpnet'.")
    method <- ifelse(is.null(object$grpnet.fit$formula), "default", "formula")
    family <- object$grpnet.fit$family$family
    if(family %in% c("multigaussian", "multinomial")){
      fit <- array(dim = c(object$grpnet.fit$nobs, length(object$grpnet.fit$ylev), 2L))
      dimnames(fit) <- list(1:object$grpnet.fit$nobs, object$grpnet.fit$ylev, c("lambda.1se", "lambda.min"))
      imp <- array(dim = c(object$grpnet.fit$ngroups - 1L, length(object$grpnet.fit$ylev), 2L))
      dimnames(imp) <- list(object$grpnet.fit$term.labels[-1], object$grpnet.fit$ylev, c("lambda.1se", "lambda.min"))
      if(method == "default"){
        fit[,,1] <- predict(object, newx = object$grpnet.fit$data$x, type = "response", s = "lambda.1se", ...)
        fit[,,2] <- predict(object, newx = object$grpnet.fit$data$x, type = "response", s = "lambda.min", ...)
        imp[,,1] <- predict(object, newx = object$grpnet.fit$data$x, type = "importance", s = "lambda.1se")
        imp[,,2] <- predict(object, newx = object$grpnet.fit$data$x, type = "importance", s = "lambda.min")
      } else {
        fit[,,1] <- predict(object, newdata = object$grpnet.fit$data, type = "response", s = "lambda.1se", ...)
        fit[,,2] <- predict(object, newdata = object$grpnet.fit$data, type = "response", s = "lambda.min", ...)
        imp[,,1] <- predict(object, newdata = object$grpnet.fit$data, type = "importance", s = "lambda.1se")
        imp[,,2] <- predict(object, newdata = object$grpnet.fit$data, type = "importance", s = "lambda.min")
      }
    } else {
      fit.type <- ifelse(family %in% c("binomail", "ordinal","svm1", "svm2"), "class", "response")
      if(method == "default"){
        fit <- cbind(predict(object, newx = object$grpnet.fit$data$x, s = "lambda.1se", type = fit.type, ...),
                     predict(object, newx = object$grpnet.fit$data$x, s = "lambda.min", type = fit.type, ...))
        imp <- cbind(predict(object, newx = object$grpnet.fit$data$x, type = "importance", s = "lambda.1se"),
                     predict(object, newx = object$grpnet.fit$data$x, type = "importance", s = "lambda.min"))
      } else {
        fit <- cbind(predict(object, newdata = object$grpnet.fit$data, s = "lambda.1se", type = fit.type, ...),
                     predict(object, newdata = object$grpnet.fit$data, s = "lambda.min", type = fit.type, ...))
        imp <- cbind(predict(object, newdata = object$grpnet.fit$data, type = "importance", s = "lambda.1se"),
                     predict(object, newdata = object$grpnet.fit$data, type = "importance", s = "lambda.min"))
      }
      colnames(fit) <- colnames(imp) <- c("lambda.1se", "lambda.min")
      if(fit.type == "class") fit <- as.data.frame(fit)
    }
    penalties <- c("LASSO", "MCP", "SCAD")
    res <- list(family = object$grpnet.fit$family$family, 
                penalty = penalties[object$grpnet.fit$args$penalty],
                nobs = object$grpnet.fit$nobs,
                ngroups = object$grpnet.fit$ngroups,
                lambda = c("lambda.1se" = object$lambda.1se, "lambda.min" = object$lambda.min),
                dev.ratio = c("lambda.1se" = object$grpnet.fit$dev.ratio[object$index[2]],
                              "lambda.min" = object$grpnet.fit$dev.ratio[object$index[1]]),
                fit = fit, 
                act = abs(imp) > 0.0,
                imp = imp)
    class(res) <- "summary.cv.grpnet"
    return(res)
    
  }  # summary.cv.grpnet

print.summary.cv.grpnet <-
  function(x, ...){
    cat("\n")
    cat("family =", x$family, "\n")
    cat("penalty =", x$penalty, "\n")
    cat("n =", x$nobs, "observations\n")
    cat("K =", x$ngroups, "groups\n")
    cat("\n% Null Deviance Explained:\n")
    print(100 * x$dev.ratio, digits = 6)
    cat("\nVariable Importance:\n")
    if(x$family %in% c("multigaussian", "multinomial")) cat("\n")
    print(100 * x$imp, digits = 6)
    if(!(x$family %in% c("multigaussian", "multinomial"))) cat("\n")
  }