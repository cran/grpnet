\name{family.grpnet}
\alias{family.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare 'family' Argument for grpnet
}
\description{
Takes in the \code{family} argument from \code{\link{grpnet}} and returns a list containing the information needed for fitting and/or tuning the model.
}
\usage{
family.grpnet(object, theta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
one of the following characters specifying the exponential family: \code{"gaussian"}, \code{"multigaussian"}, \code{"binomial"}, \code{"multinomial"}, \code{"poisson"}, \code{"negative.binomial"}, \code{"Gamma"}, \code{"inverse.gaussian"}}
\item{theta}{
Additional ("size") parameter for negative binomial responses, where the variance function is defined as \eqn{V(\mu) = \mu + \mu^2/ \theta}
}
}
\details{
There is only one available link function for each \code{family}: \cr
* gaussian (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* multigaussian (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* binomial (logit): \eqn{\log(\frac{\pi}{1 - \pi}) = \mathbf{X}^\top \boldsymbol\beta} \cr
* multinomial (symmetric): \eqn{\pi_\ell = \frac{\exp(\mathbf{X}^\top \boldsymbol\beta_\ell)}{\sum_{l = 1}^m \exp(\mathbf{X}^\top \boldsymbol\beta_l)}} \cr
* poisson (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* negative.binomial (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* Gamma (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* inverse.gaussian (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
}
\value{
List with components:
\item{family }{same as input object, i.e., character specifying the family}
\item{linkinv }{function for computing inverse of link function}
\item{dev.resids }{function for computing deviance residuals}
}
\references{
Helwig, N. E. (2024). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics}. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For \code{gaussian} family, this returns the full output produced by \code{\link[stats]{gaussian}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{grpnet}} for fitting group elastic net regularization paths

\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}
}
\examples{
family.grpnet("gaussian")

family.grpnet("multigaussian")

family.grpnet("binomial")

family.grpnet("multinomial")

family.grpnet("poisson")

family.grpnet("negbin", theta = 10)

family.grpnet("Gamma")

family.grpnet("inverse.gaussian")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ models }
\keyword{ regression }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
