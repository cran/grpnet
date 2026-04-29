summary.grpnet <-
  function(object, ...){
    # summarize grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2026-04-28
    
    if(!inherits(object, "grpnet")) stop("Input 'object' should be of class 'grpnet'.")
    method <- ifelse(is.null(object$formula), "default", "formula")
    if(method == "default"){
      fit <- predict(object, newx = object$data$x, ...)
      imp <- predict(object, newx = object$data$x, type = "importance")
    } else {
      fit <- predict(object, newdata = object$data, ...)
      imp <- predict(object, newdata = object$data, type = "importance")
    }
    if(object$family$family %in% c("multigaussian", "multinomial")){
      imp0 <- imp
      imp <- array(dim = c(nrow(imp0[[1]]), length(imp0), ncol(imp0[[1]])),
                   dimnames = list(rownames(imp0[[1]]), names(imp0), colnames(imp0[[1]])))
      for(k in 1:length(imp0)) imp[,k,] <- imp0[[k]]
    }
    penalties <- c("LASSO", "MCP", "SCAD")
    res <- list(family = object$family$family, 
                penalty = penalties[object$args$penalty],
                nobs = object$nobs,
                ngroups = object$ngroups,
                lambda = object$lambda,
                dev.ratio = object$dev.ratio,
                fit = fit, 
                act = abs(imp) > 0.0,
                imp = imp)
    class(res) <- "summary.grpnet"
    return(res)
  } # summary.grpnet

print.summary.grpnet <-
  function(x, ...){
    cat("\n")
    cat("family =", x$family, "\n")
    cat("penalty =", x$penalty, "\n")
    cat("n =", x$nobs, "observations\n")
    cat("K =", x$ngroups, "groups\n")
    nlam <- length(x$lambda)
    lamseq <- seq(1, nlam, length.out = min(5, nlam))
    cat("\n% Null Deviance Explained:\n")
    names(x$dev.ratio) <- paste0("s", 1:length(x$lambda))
    dr <- x$dev.ratio[lamseq]
    print(100 * dr, digits = 6)
    cat("\nVariable Importance:\n")
    if(x$family %in% c("multigaussian", "multinomial")){
      respvars <- dimnames(x$imp)[[2]]
      cat("\n")
      for(k in 1:length(respvars)){
        cat("  ", respvars[k], ":\n", sep = "")
        print(100 * x$imp[,k,lamseq], digits = 6)
        cat("\n")
      }
    } else {
      print(100 * x$imp[,lamseq], digits = 6)
    }
    cat("\n")
  }