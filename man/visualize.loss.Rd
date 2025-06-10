\name{visualize.loss}
\alias{visualize.loss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots grpnet Loss Function on Link or Response Scale
}
\description{
Makes a plot or returns a data frame containing the specified loss function evaluated at a sequence of input values.
}
\usage{
visualize.loss(x = seq(-3, 3, length.out = 1001),
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
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  sequence of linear predictor values at which to evaluate the loss.
  }
  \item{family}{
  Character specifying the assumed distribution for the response variable. Partial matching is allowed. See options below.
  }
  \item{theta}{
  For SVM1: additional ("smoothing") parameter, that controls the smoothing rate of the hinge loss function.
For negative binomial: additional ("size") parameter, where the variance function is defined as \eqn{V(\mu) = \mu + \mu^2/ \theta}
  }
  \item{type}{
  Default of \code{type = "link"} plots \code{x} versus the loss. Change to \code{type = "response"} to plot \code{mu} versus the loss.
  }
  \item{y}{
  Response value used to compute loss. Note that the loss function is interpreted as a function of \code{x} given \code{y}.
  }
  \item{plot}{
  if \code{TRUE} (default), then the result is plotted; otherwise the result is returned as a data frame.
  }
  \item{add}{
  if \code{FALSE} (default), then the \code{\link{plot}} function is used to create a new graphic; otherwise the \code{\link{lines}} function is used to add lines to an existing graphic.
  }
  \item{...}{
  additional arguments passed to \code{plot} function, e.g., \code{xlim}, \code{ylim}, etc.
  }
}
\details{
\bold{grpnet} implements the following loss functions:
\describe{
\item{gaussian/multigaussian}{\eqn{L = (y - \mu)^2}}
\item{svm1}{\eqn{L = \left\{ \begin{array}{ll}
\frac{1}{2\theta}(1 - \mu y)_{+}^2  & 1 - \theta < \mu y \\
1 - \mu y - \theta / 2 & \mu y \leq 1 - \theta 
\end{array} \right.}}
\item{svm2}{\eqn{L = (1 - \mu y)_{+}^2}}
\item{logit}{\eqn{L = \log(1 + \exp(-\mu y))}}
\item{binomial}{\eqn{L = -y\log(\mu) - (1-y) \log(1-\mu)}}
\item{multinomial}{\eqn{L = -\sum_{l=1}^m I(y=l) \log(\mu_l)}}
\item{poisson}{\eqn{L = \mu - y \log(\mu)}}
\item{negative.binomial}{L = \eqn{(\theta + y)  \log(\theta + \mu) - y  \log(\mu) + c} \cr where \eqn{c = \log(\Gamma(\theta)) - \log(\Gamma(\theta + y)) - \theta  \log(\theta)} \cr is a constant with respect to \eqn{\mu}}
\item{Gamma}{\eqn{L = \log(\mu) + y / \mu}}
\item{inverse.gaussian}{\eqn{L = (y - \mu)^2 / (\mu^2 y)}}
}
}

\value{
If \code{plot = TRUE}, a plot is produced and nothing is returned.

If \code{plot = FALSE}, a data frame is returned with columns:
\item{eta}{linear predictor, i.e., fitted values on link scale (same as input \code{x})}
\item{mu}{expected value, i.e., fitted values on response scale (\code{family$linkinv(x)})}
\item{loss}{loss function evaluation for given \code{eta} (and \code{y})}
}
\references{
Helwig, N. E. (2025). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics, 34}(1), 239-252. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{visualize.penalty}} for plotting penalty function

\code{\link{visualize.shrink}} for plotting shrinkage operator
}
\examples{
visualize.loss(family = "gaussian")

visualize.loss(family = "svm1", theta = 0.1)

visualize.loss(family = "svm2")

visualize.loss(family = "logit")

visualize.loss(family = "binomial")

visualize.loss(family = "poisson")

visualize.loss(family = "negative.binomial", theta = 10)

visualize.loss(family = "Gamma")

visualize.loss(family = "inverse.gaussian")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ models }
\keyword{ regression }
\keyword{ dplot }
\keyword{ hplot }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
