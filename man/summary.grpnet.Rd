\name{summary.grpnet}
\alias{summary.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize grpnet Fits
}
\description{
Obtain fitted values, variable activations, variable importances from a fit group elastic net (grpnet) object.
}
\usage{
\method{summary}{grpnet}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "grpnet"
}
  \item{...}{
Additional arguments passed to \code{\link{predict.grpnet}} function.
}
}
\details{
Passes the training data to the \code{\link{predict.grpnet}} function twice: (1) to calculate fitted values, and (2) to calculate variable importances.
}
\value{
\item{family}{Name of exponential family/loss function.}
\item{penalty}{Name of L1 group penalty.}
\item{nobs}{Number of observations.}
\item{ngroups}{Number of coefficient groups.}
\item{lambda}{Vector of lambda values.}
\item{dev.ratio}{Vector of proportions of null deviance explained.}
\item{fit}{Fitted values obtained from \code{\link{predict.grpnet}}.}
\item{act}{Variable activations obtained from \code{\link{predict.grpnet}}.}
\item{imp}{Variable importances obtained from \code{\link{predict.grpnet}}.}
}
\references{
Helwig, N. E. (2025). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics, 34}(1), 239-252. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The variables activations are defined as:  \code{act <- abs(imp) > 0.0}

For most response families, the returned components \code{fit}, \code{act}, and \code{imp} are matrices of dimension \code{K} by \code{nlambda}.

For families with multivariate responses (i.e., \code{multigaussian} and \code{multinomial}), these components are arrays of dimension \code{K} by \code{L} by \code{nlambda} where \code{L} is the dimension of the response.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{grpnet}} for fitting grpnet regularization paths
}
\examples{
# see 'grpnet' for summary examples
?grpnet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
