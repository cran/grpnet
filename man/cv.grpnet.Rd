\name{cv.grpnet}
\alias{cv.grpnet}
\alias{cv.grpnet.default}
\alias{cv.grpnet.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation for grpnet
}
\description{
Implements k-fold cross-validation for \code{\link{grpnet}} to find the regularization parameters that minimize the prediction error (deviance, mean squared error, mean absolute error, or misclassification rate).
}
\usage{
cv.grpnet(x, \dots)

\method{cv.grpnet}{default}(x, 
          y, 
          group,
          weights = NULL,
          offset = NULL,
          alpha = c(0.01, 0.25, 0.5, 0.75, 1),
          gamma = c(3, 4, 5),
          type.measure = NULL,
          nfolds = 10, 
          foldid = NULL,
          same.lambda = FALSE,
          parallel = FALSE, 
          cluster = NULL, 
          verbose = interactive(), 
          adaptive = FALSE,
          power = 1,
          ...)
           
\method{cv.grpnet}{formula}(formula,
          data, 
          use.rk = TRUE,
          weights = NULL,
          offset = NULL,
          alpha = c(0.01, 0.25, 0.5, 0.75, 1),
          gamma = c(3, 4, 5),
          type.measure = NULL,
          nfolds = 10, 
          foldid = NULL, 
          same.lambda = FALSE,
          parallel = FALSE, 
          cluster = NULL, 
          verbose = interactive(), 
          adaptive = FALSE,
          power = 1,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Model (design) matrix of dimension \code{nobs} by \code{nvars} (\eqn{n \times p}).
}
  \item{y}{
Response vector of length \eqn{n} or matrix of dimension \eqn{n \times m}. Note that matrix inputs are (i) required for multigaussian family, (ii) allowed for binomial and multinomial families (see "Binomial and multinomial" section in \code{\link{grpnet}}), and (iii) not permitted for other families.
}
  \item{group}{
Group label vector (factor, character, or integer) of length \eqn{p}. Predictors with the same label are grouped together for regularization.
}
  \item{formula}{
Model formula: a symbolic description of the model to be fitted. Uses the same syntax as \code{\link{lm}} and \code{\link{glm}}. 
}
  \item{data}{
Optional data frame containing the variables referenced in \code{formula}.
}
  \item{use.rk}{
If \code{TRUE} (default), the \code{\link{rk.model.matrix}} function is used to build the model matrix. Otherwise, the \code{\link{model.matrix}} function is used to build the model matrix. Additional arguments to the \code{\link{rk.model.matrix}} function can be passed via the \code{...} argument.
}
  \item{weights}{
Optional vector of length \eqn{n} with non-negative weights to use for weighted (penalized) likelihood estimation. Defaults to a vector of ones.
}
  \item{offset}{
Optional vector of length \eqn{n} with an a priori known term to be included in the model's linear predictor. Defaults to a vector of zeros.  
}
  \item{alpha}{
Scalar or vector specifying the elastic net tuning parameter \eqn{\alpha}. If \code{alpha} is a vector (default), then (a) the same \code{foldid} is used to compute the cross-validation error for each \eqn{\alpha}, and (b) the solution for the optimal \eqn{\alpha} is returned. 
}
\item{gamma}{
Scalar or vector specifying the penalty hyperparameter \eqn{\gamma} for MCP or SCAD. If \code{gamma} is a vector (default), then (a) the same \code{foldid} is used to compute the cross-validation error for each \eqn{\gamma}, and (b) the solution for the optimal \eqn{\gamma} is returned. 
}   
  \item{type.measure}{
Loss function for cross-validation. Options include: \code{"deviance"} for model deviance, \code{"mse"} for mean squared error, \code{"mae"} for mean absolute error, or \code{"class"} for classification error. Note that \code{"class"} is only available for binomial and multinomial families. The default is classification error (for binomial and multinomial) or mean absolute error (others).
}
  \item{nfolds}{
Number of folds for cross-validation. 
}
  \item{foldid}{
Optional vector of length \eqn{n} giving the fold identification for each observation. Must be coercible into a factor. After coersion, the \code{nfolds} argument is defined as \code{nfolds = nlevels(foldid)}.
}
  \item{same.lambda}{
Logical specfying if the same \eqn{\lambda} sequence should be used for fitting the model to each fold's data. If \code{FALSE} (default), the \eqn{\lambda} sequence is determined separately holding out each fold, and the \eqn{\lambda} sequence from the full model is used to align the predictions. If \code{TRUE}, the \eqn{\lambda} sequence from the full model is used to fit the model for each fold. The default often provides better (i.e., more stable) computational performance.
}
  \item{parallel}{
Logical specifying if sequential computing (default) or parallel computing should be used. If \code{TRUE}, the fitting for each fold is parallelized.
}
  \item{cluster}{
Optional cluster to use for parallel computing. If \code{parallel = TRUE} and \code{cluster = NULL}, then the cluster is defined \code{cluster = makeCluster(2L)}, which uses two cores. Recommended usage: \code{cluster = makeCluster(detectCores())}
}
  \item{verbose}{
Logical indicating if the fitting progress should be printed. Defaults to \code{TRUE} in interactive sessions and \code{FALSE} otherwise. 
}
\item{adaptive}{
Logical indicating if the adaptive group elastic net should be used (see Note).
}
\item{power}{
If \code{adaptive = TRUE}, then the adaptive penalty weights are defined by dividing the original penalty weights by \code{tapply(coef, group, norm, type = "F")^power}.
}
  \item{...}{
Optional additional arguments for \code{\link{grpnet}} (e.g., \code{standardize}, \code{penalty.factor}, etc.)
}
}
\details{
This function calls the \code{\link{grpnet}} function \code{nfolds+1} times: once on the full dataset to obtain the \code{lambda} sequence, and once holding out each fold's data to evaluate the prediction error. The syntax of (the default S3 method for) this function closely mimics that of the \code{cv.glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010). 

Let \eqn{\mathbf{D}_u = \{\mathbf{y}_u, \mathbf{X}_u\}} denote the \eqn{u}-th fold's data, let \eqn{\mathbf{D}_{[u]} = \{\mathbf{y}_{[u]}, \mathbf{X}_{[u]}\}} denote the full dataset excluding the \eqn{u}-th fold's data, and let \eqn{\boldsymbol\beta_{\lambda [u]}} denote the coefficient estimates obtained from fitting the model to \eqn{\mathbf{D}_{[u]}} using the regularization parameter \eqn{\lambda}. 

The cross-validation error for the \eqn{u}-th fold is defined as
\deqn{E_u(\lambda) = C(\boldsymbol\beta_{\lambda [u]} , \mathbf{D}_u)}
where \eqn{C(\cdot , \cdot)} denotes the cross-validation loss function that is specified by \code{type.measure}. For example, the \code{"mse"} loss function is defined as
\deqn{C(\boldsymbol\beta_{\lambda [u]} , \mathbf{D}_u) = \| \mathbf{y}_u - \mathbf{X}_u \boldsymbol\beta_{\lambda [u]} \|^2}
where \eqn{\| \cdot \|} denotes the L2 norm.

The mean cross-validation error \code{cvm} is defined as
\deqn{\bar{E}(\lambda) = \frac{1}{v} \sum_{u = 1}^v E_u(\lambda) }
where \eqn{v} is the total number of folds. The standard error \code{cvsd} is defined as
\deqn{S(\lambda) = \sqrt{ \frac{1}{v (v - 1)} \sum_{u=1}^v (E_u(\lambda) - \bar{E}(\lambda))^2 } }
which is the classic definition of the standard error of the mean. 
}
\value{
\item{lambda}{regularization parameter sequence for the full data}
\item{cvm}{mean cross-validation error for each \code{lambda}}
\item{cvsd}{estimated standard error of \code{cvm}}
\item{cvup}{upper curve: \code{cvm + cvsd}}
\item{cvlo}{lower curve: \code{cvm - cvsd}}
\item{nzero}{number of non-zero groups for each \code{lambda}}
\item{grpnet.fit}{fitted grpnet object for the full data}
\item{lambda.min}{value of \code{lambda} that minimizes \code{cvm}}
\item{lambda.1se}{largest \code{lambda} such that \code{cvm} is within one \code{cvsd} from the minimum (see Note)}
\item{index}{two-element vector giving the indices of \code{lambda.min} and \code{lambda.1se} in the \code{lambda} vector, i.e., \code{c(minid, se1id)} as defined in the Note}
\item{type.measure}{loss function for cross-validation (used for plot label)}
\item{call}{matched call}
\item{time}{runtime in seconds to perform k-fold CV tuning}
\item{tune}{data frame containing the tuning results, i.e., min(cvm) for each combination of \code{alpha} and/or \code{gamma}}
}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}

Helwig, N. E. (2025). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics, 34}(1), 239-252. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
When \code{adaptive = TRUE}, the adaptive group elastic net is used: \cr
(1) an initial fit with \code{alpha = 0} estimates the \code{penalty.factor} \cr
(2) a second fit using estimated \code{penalty.factor} is returned \cr

\code{lambda.1se} is defined as follows: \cr
\code{minid <- which.min(cvm)} \cr
\code{min1se <- cvm[minid] + cvsd[minid]} \cr
\code{se1id <- which(cvm <= min1se)[1]} \cr
\code{lambda.1se <- lambda[se1id]}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.cv.grpnet}} for plotting the cross-validation error curve

\code{\link{predict.cv.grpnet}} for predicting from \code{\link{cv.grpnet}} objects

\code{\link{grpnet}} for fitting group elastic net regularization paths
}
\examples{
\donttest{
######***######   family = "gaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = mpg)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto)

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "multigaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = (mpg, displacement))
y <- as.matrix(auto[,c(1,3)])
set.seed(1)
mod <- cv.grpnet(y ~ ., data = auto[,-c(1,3)], family = "multigaussian",
                 standardize.response = TRUE)

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "binomial"   ######***######

# load data
data(auto)

# redefine origin (Domestic vs Foreign)
auto$origin <- ifelse(auto$origin == "American", "Domestic", "Foreign")

# 10-fold cv (default method, response = origin with 2 levels)
set.seed(1)
mod <- cv.grpnet(origin ~ ., data = auto, family = "binomial")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "multinomial"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin with 3 levels)
set.seed(1)
mod <- cv.grpnet(origin ~ ., data = auto, family = "multinomial")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "poisson"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = horsepower)
set.seed(1)
mod <- cv.grpnet(horsepower ~ ., data = auto, family = "poisson")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "negative.binomial"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = horsepower)
set.seed(1)
mod <- cv.grpnet(horsepower ~ ., data = auto, family = "negative.binomial")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "Gamma"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto, family = "Gamma")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)



######***######   family = "inverse.gaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto, family = "inverse.gaussian")

# print min and 1se solution info
mod

# plot cv error curve
plot(mod)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
