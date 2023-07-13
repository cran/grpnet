# Some code (for printing and plotting) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

grpnet <- 
  function(x, ...){
    UseMethod("grpnet")
  } # end grpnet

print.grpnet <-
  function(x, ...){
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    df <- data.frame(x$nzgrp,
                     x$df, 
                     round(100 * x$dev.ratio, 2),
                     x$lambda)
    colnames(df) <- c("nGrp", "Df", "%Dev", "Lambda")
    intercept <- ifelse(max(abs(x$a0)) > 0, TRUE, FALSE)
    print(df, digits = 6)
    cat("\n")
  } # end print.grpnet

plot.grpnet <-
  function(x, color.by.group = TRUE, col = NULL, ...){
    if(color.by.group) {
      if(is.null(col)){
        col <- 1:x$ngroups
      } else {
        if(length(col) != x$ngroups) stop("Input 'col' must be of length x$ngroups")
      }
      colors <- col[as.integer(as.factor(x$group))]
    } else {
      if(is.null(col)) col <- "black"
      colors <- rep(col[1], nrow(x$beta))
    }
    if(x$family$family == "multinomial"){
      for(j in 1:length(x$beta)){
        plot(log(x$lambda), x$beta[[j]][1,], ylim = extendrange(x$beta),
             xlab = "Log Lambda", ylab = "Coefficients", t = "n", ...)
        legend("top", legend = x$ylev[j], bty = "n", cex = 0.8)
        for(k in 1:nrow(x$beta[[j]])) {
          lines(log(x$lambda), x$beta[[j]][k,], col = colors[k])
        }
      }
    } else {
      plot(log(x$lambda), x$beta[1,], ylim = extendrange(x$beta),
           xlab = "Log Lambda", ylab = "Coefficients", t = "n", ...)
      for(k in 1:nrow(x$beta)) {
        lines(log(x$lambda), x$beta[k,], col = colors[k])
      }
    }
    abline(h = 0)
  } # end plot.grpnet