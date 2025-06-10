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
  two options: (1) an object of class "grpnet" or "cv.grpnet"; or (2) a character specifying the exponential family: \code{"gaussian"}, \code{"multigaussian"}, \code{"svm1"}, \code{"svm2"}, \code{"logit"}, \code{"binomial"}, \code{"multinomial"}, \code{"poisson"}, \code{"negative.binomial"}, \code{"Gamma"}, \code{"inverse.gaussian"}}
\item{theta}{
positive scalar that serves as an additional hyperparameter for various loss functions.

svm1: additional parameter that controls the smoothing rate for the hinge loss function (see Note below).

negative.binomial: size parameter such that the variance function is defined as \eqn{V(\mu) = \mu + \mu^2/ \theta}
}
}
\details{
There is only one available link function for each \code{family}: \cr
* gaussian (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* multigaussian (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* svm1/svm2 (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* binomial/logit (logit): \eqn{\log(\frac{\pi}{1 - \pi}) = \mathbf{X}^\top \boldsymbol\beta} \cr
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
Helwig, N. E. (2025). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics, 34}(1), 239-252. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For \code{gaussian} family, this returns the full output produced by \code{\link[stats]{gaussian}}.

For \code{svm1} family, the quadratically smoothed hinge loss is defined as
\deqn{
\mathrm{svm1}(z) = \left\{ \begin{array}{ll}
0 & z > 1 \\
(1 - z)^2 / (2 \theta) & 1 - \theta < z \leq 1 \\
1 - z - \theta / 2 & z \leq 1 - \theta \\
\end{array} \right.
}
where \eqn{z = Y \eta} with \eqn{Y \in \{-1,1\}} denoting the response and \eqn{\eta = \mathbf{X}^\top \boldsymbol\beta} denoting the linear predictor. Note that the \code{svm1} loss function approaches the support vector machine (i.e., hinge) loss function as \eqn{\theta \rightarrow 0}.

For \code{svm2} family, the squared hinge loss is defined as
\deqn{
\mathrm{svm2}(z) = \left\{ \begin{array}{ll}
0 & z > 1 \\
(1 - z)^2 & z \leq 1 \\
\end{array} \right.
}
where \eqn{z = Y \eta} with \eqn{Y \in \{-1,1\}} denoting the response and \eqn{\eta = \mathbf{X}^\top \boldsymbol\beta} denoting the linear predictor. Note that the \code{svm1} loss function approaches the support vector machine (i.e., hinge) loss function as \eqn{\theta \rightarrow 0}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{visualize.loss}} for plotting loss functions

\code{\link{grpnet}} for fitting group elastic net regularization paths

\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}
}
\examples{
family.grpnet("gaussian")

family.grpnet("multigaussian")

family.grpnet("svm1", theta = 0.1)

family.grpnet("svm2")

family.grpnet("logit")

family.grpnet("binomial")

family.grpnet("multinomial")

family.grpnet("poisson")

family.grpnet("negative.binomial", theta = 10)

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
