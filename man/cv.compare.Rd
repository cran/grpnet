\name{cv.compare}
\alias{cv.compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare Multiple cv.grpnet Solutions
}
\description{
Creates a plot (default) or returns a data frame (otherwise) that compares the cross-validation error for multiple \code{\link{cv.grpnet}} fits.
}
\usage{
cv.compare(x,
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
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a single \code{\link{cv.grpnet}} object or a list of \code{\link{cv.grpnet}} objects.
}
  \item{s}{
the tuning parameter value at which to plot results (if \code{x} is a list).
}
  \item{plot}{
switch controlling whether a plot is produced (default) versus data frame.
}
  \item{at}{
x-axis coordinates for plotting the cv error for each solution.
}
  \item{nse}{
number of standard errors to use for error bars in plot.
}
  \item{point.col}{
color for point used to plot the average of the cv error.
}
  \item{line.col}{
color for lines used to plot the standard error for the cv error.
}
  \item{lwd}{
width of lines used to plot the standard error for the cv error.  
}
  \item{bwd}{
width of standard error bars in terms of proportion of \code{range(x)}.
}
  \item{labels}{
labels for x-axis tick marks. Defaults to \code{names(x)}.
}
  \item{xlim}{
axis limits for abscissa (x-axis)
}
  \item{ylim}{
axis limits for ordinate (y-axis)
}
  \item{xlab}{
axis label for abscissa (x-axis)
}
  \item{ylab}{
axis label for ordinate (y-axis)
}
  \item{...}{
additional arguments passed to plotting functions.  
}

}
\details{
Default behavior creates a plot that displays the mean cv error +/- 1 se for each of the requested solutions. 

If the input \code{x} is a single \code{\link{cv.grpnet}} object, then the function plots the lambda.min and lambda.1se solutions. 

If the input \code{x} is a list of \code{\link{cv.grpnet}} objects, then the function plots either the lambda.min or the lambda.1se solution (controlled by \code{s} argument) for all of the input models.
}
\value{
When \code{plot = TRUE}, there is no return value (it produces a plot)

When \code{plot = FALSE}, a data.frame is returned with the mean cv error (and se) for each solution
}
\references{
Helwig, N. E. (2024). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics}. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.cv.grpnet}} for plotting cv error path (for all lambdas)

\code{\link{plot.grpnet}} for plotting regularization path (for single lambda)
}
\examples{
# load data
data(auto)

# LASSO penalty
set.seed(1)
mod1 <- cv.grpnet(mpg ~ ., data = auto, nfolds = 5, alpha = 1)

# MCP penalty
set.seed(1)
mod2 <- cv.grpnet(mpg ~ ., data = auto, nfolds = 5, alpha = 1, penaly = "MCP")

# SCAD penalty
set.seed(1)
mod3 <- cv.grpnet(mpg ~ ., data = auto, nfolds = 5, alpha = 1, penaly = "SCAD")

# compare lambda.min and lambda.1se for mod1
cv.compare(mod1)

# compare lambda.1se for mod1, mod2, mod3
cv.compare(x = list(mod1, mod2, mod3), labels = c("LASSO", "MCP", "SCAD"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ graphs }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
