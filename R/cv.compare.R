cv.compare <-
  function(x,
           s = c("lambda.1se", "lambda.min"), 
           plot = TRUE, 
           at = 1:length(x),
           nse = 1, 
           point.col = "red", 
           line.col = "gray", 
           lwd = 2, 
           bwd = 0.02,
           labels = NULL,
           xlim = NULL,
           ylim = NULL,
           xlab = NULL,
           ylab = NULL,
           ...){
    # compare multiple cv.grpnet prediction errors
    # Nathaniel E. Helwig (helwig@umn.edu)
    # 2024-06-05
    
    
    ### xlist
    x <- as.list(x)
    nx <- length(x)
    
    ### check x class
    xcheck <- inherits(x, "cv.grpnet")
    
    ### is x cv.grpnet object?
    if(xcheck){
      
      xnames <- c("lambda.min", "lambda.1se")
      at <- 1:2
      rx <- 1
      cvm <- x$cvm[x$index]
      cvsd <- x$cvsd[x$index]
      type.measure <- x$type.measure
      family <- x$grpnet.fit$family$family
      
    } else {
      
      ### check x
      xcheck <- sapply(x, inherits, what = "cv.grpnet")
      if(any(!xcheck)) stop("Input 'x' must be a list of 'cv.grpnet' objects")
      
      ## xnames
      xnames <- labels
      if(is.null(xnames)) xnames <- paste0("m", 1:nx)
      
      ## check s
      s <- as.character(s[1])
      s <- pmatch(s, c("lambda.min", "lambda.1se"))
      if(is.na(s)) stop("Input 's' must be 'lambda.min' or 'lambda.1se'")
      s <- c("lambda.min", "lambda.1se")[s]
      sid <- ifelse(s == "lambda.min", 1, 2)
      
      ## check at
      at <- as.numeric(at)
      if(length(at) != nx) stop("Input 'at' must satisfy: length(x) == length(at)")
      rx <- range(at)
      rx <- rx[2] - rx[1]
      
      ## extract mean cv error and sd
      cvm <- cvsd <- rep(NA, nx)
      for(i in 1:nx){
        index <- x[[i]]$index[sid]
        cvm[i] <- x[[i]]$cvm[index]
        cvsd[i] <- x[[i]]$cvsd[index]
      }
      
      ## get type.measure andfamily
      type.measure <- x[[1]]$type.measure
      family <- x[[1]]$grpnet.fit$family$family
      
    } # end if(xcheck)
    
    ### return data.frame?
    if(!plot) return(data.frame(model = xnames, cvm = cvm, cvsd = cvsd))
    
    ### check xlab and ylab
    if(is.null(xlab)) xlab <- "Models"
    if(is.null(ylab)){
      if(type.measure == "deviance"){
        titlecase <- function(x) {
          cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), 
                                     tolower(substr(x, 2, nchar(x))))
          words <- strsplit(x, split = " ", fixed = TRUE)
          words <- sapply(words, cap1)
          paste(words, collapse = " ")
        }
        ylab <- paste0(titlecase(family), "Deviance")
      } else if(type.measure == "mse"){
        ylab <- "Mean Squared Error"
      } else if(type.measure == "mae"){
        ylab <- "Mean Absolute Error"
      } else if(type.measure == "class"){
        ylab <- "Misclassification Error"
      }
    }
    
    ### check xlim ylim
    if(is.null(xlim)){
      xlim <- extendrange(at)
      xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
    }
    if(is.null(ylim)){
      ylim <- extendrange(c(cvm - nse*cvsd, cvm + nse*cvsd))
    }
    
    ### make plot
    plot(at, cvm, pch = 19, axes = F, xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, col = point.col, type = "n", ...)
    segments(x0 = at, y0 = cvm - nse*cvsd, y1 = cvm + nse*cvsd, lwd = lwd, col = line.col)
    segments(x0 = at - bwd*rx, y0 = cvm - nse*cvsd, x1 = at + bwd*rx, lwd = lwd, col = line.col)
    segments(x0 = at - bwd*rx, y0 = cvm + nse*cvsd, x1 = at + bwd*rx, lwd = lwd, col = line.col)
    points(at, cvm, pch = 19, col = point.col)
    axis(1, at = at, labels = xnames, ...)
    axis(2, ...)
    box()
    
  } # end cv.compare