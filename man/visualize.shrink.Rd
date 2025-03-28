\name{visualize.shrink}
\alias{visualize.shrink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots grpnet Shrinkage Operator or its Estimator
}
\description{
Makes a plot or returns a data frame containing the group elastic net shrinkage operator (or its estimator) evaluated at a sequence of input values.
}
\usage{
visualize.shrink(x = seq(-5, 5, length.out = 1001), 
                penalty = c("LASSO", "MCP", "SCAD"), 
                alpha = 1, 
                lambda = 1, 
                gamma = 4, 
                fitted = FALSE,
                plot = TRUE,
                subtitle = TRUE,
                legend = TRUE,
                location = "top",
                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  sequence of values at which to evaluate the penalty.
  }
  \item{penalty}{
  which penalty or penalties should be plotted?  
  }
  \item{alpha}{
  elastic net tuning parameter (between 0 and 1).  
  }
  \item{lambda}{
  overall tuning parameter (non-negative).
  }
  \item{gamma}{
  additional hyperparameter for MCP (>1) or SCAD (>2).
  }
  \item{fitted}{
  if \code{FALSE} (default), then the shrinkage operator is plotted; otherwise the shrunken estimator is plotted.
  }
  \item{plot}{
  if \code{TRUE} (default), then the result is plotted; otherwise the result is returned as a data frame.
  }
  \item{subtitle}{
  if \code{TRUE} (default), then the hyperparameter values are displayed in the subtitle. 
  }
  \item{legend}{
  if \code{TRUE} (default), then a legend is included to distinguish the different \code{penalty} types.
  }
  \item{location}{
  the legend's location; ignored if \code{legend = FALSE}.
  }
  \item{...}{
  addition arguments passed to \code{plot} function, e.g., \code{xlim}, \code{ylim}, etc.
  }
}
\details{
The updates for the group elastic net estimator have the form
\deqn{\boldsymbol\beta_{\alpha, \lambda}^{(t+1)} = S_{\lambda_1, \lambda_2}(\|\mathbf{b}_{\alpha, \lambda}^{(t+1)}\|) \mathbf{b}_{\alpha, \lambda}^{(t+1)}}
where \eqn{S_{\lambda_1, \lambda_2}(\cdot)} is a shrinkage and selection operator, and 
\deqn{\mathbf{b}_{\alpha, \lambda}^{(t+1)} = \boldsymbol\beta_{\alpha, \lambda}^{(t)} + (\delta_{(t)} \epsilon)^{-1} \mathbf{g}^{(t)}}
is the unpenalized update with \eqn{\mathbf{g}^{(t)}} denoting the current gradient.


Note that \eqn{\lambda_1 = \lambda \alpha} is the L1 tuning parameter, \eqn{\lambda_2 = \lambda (1-\alpha)} is the L2 tuning parameter, \eqn{\delta_{(t)}} is an upper-bound on the weights appearing in the Fisher information matrix, and \eqn{\epsilon} is the largest eigenvalue of the Gramm matrix \eqn{n^{-1} \mathbf{X}^\top \mathbf{X}}. 
}
\value{
If \code{plot = TRUE}, then produces a plot.

If \code{plot = FALSE}, then returns a data frame.
}
\references{
Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. \emph{Journal of the American Statistical Association, 96}(456), 1348-1360. \doi{10.1198/016214501753382273}

Helwig, N. E. (2025). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics, 34}(1), 239-252. \doi{10.1080/10618600.2024.2362232}

Tibshirani, R. (1996). Regression and shrinkage via the Lasso. \emph{Journal of the Royal Statistical Society, Series B, 58}, 267-288. \doi{10.1111/j.2517-6161.1996.tb02080.x}

Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty. \emph{The Annals of Statistics, 38}(2), 894-942. \doi{10.1214/09-AOS729}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{visualize.penalty}} for plotting penalty function
}
\examples{
# plot shrinkage operator
visualize.shrink()

# plot shrunken estimator
visualize.shrink(fitted = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ dplot }
\keyword{ hplot }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
