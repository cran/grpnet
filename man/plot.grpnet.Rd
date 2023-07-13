\name{plot.grpnet}
\alias{plot.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Coefficients for grpnet Fits
}
\description{
Creates a profile plot of the coefficient paths for a fit group elastic net regularized GLM (grpnet) object.
}
\usage{
\method{plot}{grpnet}(x, color.by.group = TRUE, col = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "grpnet"
}
  \item{color.by.group}{
If \code{TRUE} (default), the coefficient paths are colored according to their group membership using the colors in \code{col}. If \code{FALSE}, all coefficient paths are plotted the same color.
}
  \item{col}{
If \code{color.by.group = TRUE}, this should be a vector of length \eqn{K} giving a color label for each group. If \code{color.by.group = FASLE}, this should be a character specifying a single (common) color. Default of \code{col = NULL} is the same as \code{col = 1:K} or \code{col = "black"}.
}
  \item{...}{
Additional arguments passed to the \code{\link{plot}} function.
}
}
\details{
Syntax and functionality were modeled after the \code{plot.glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010).
}
\value{
Produces a profile plot showing the value of each coefficient (y-axis) as a function of \code{log(lambda)} (x-axis). 
}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If \code{x} is a multinomial model, the coefficients for each response class are plotted in a separate plot.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{grpnet}} for fitting grpnet regularization paths

\code{\link{plot.cv.grpnet}} for plotting \code{\link{cv.grpnet}} objects
}
\examples{
# see 'grpnet' for plotting examples
?grpnet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ graphs }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
