\name{predict.cv.grpnet}
\alias{predict.cv.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Method for cv.grpnet Fits
}
\description{
Obtain predictions from a cross-validated group elastic net regularized GLM (cv.grpnet) object.
}
\usage{
\method{predict}{cv.grpnet}(object, 
        newx,
        newdata,
        s = c("lambda.1se", "lambda.min"),
        type = c("link", "response", "class", "terms", 
                 "importance", "coefficients", "nonzero", "groups", 
                 "ncoefs", "ngroups", "norm", "znorm"),
        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "cv.grpnet"
}
  \item{newx}{
Matrix of new \code{x} scores for prediction (default S3 method). Must have \eqn{p} columns arranged in the same order as the \code{x} matrix used to fit the model.
}
  \item{newdata}{
Data frame of new \code{data} scores for prediction (S3 "formula" method). Must contain all variables in the \code{formula} used to fit the model.
}
  \item{s}{
Lambda value(s) at which predictions should be obtained. Can input a character ("lambda.min" or "lambda.1se") or a numeric vector. Default of "lambda.min" uses the \code{lambda} value that minimizes the mean cross-validated error. 
}
  \item{type}{
Type of prediction to return. "link" gives predictions on the link scale (\eqn{\eta}). "response" gives predictions on the mean scale (\eqn{\mu}). "class" gives predicted class labels (for "binomial" and "multinomial" families). "terms" gives the predictions for each term (group) in the model (\eqn{\eta_k}). "importance" gives the variable importance index for each term (group) in the model. "coefficients" returns the coefficients used for predictions. "nonzero" returns a list giving the indices of non-zero coefficients for each \code{s}. "groups" returns a list giving the labels of non-zero groups for each \code{s}. "ncoefs" returns the number of non-zero coefficients for each \code{s}. "ngroups" returns the number of non-zero groups for each \code{s}. "norm" returns the L2 norm of each group's (raw) coefficients for each \code{s}. "znorm" returns the L2 norm of each group's standardized coefficients for each \code{s}.
}
  \item{...}{
Additional arguments (ignored)
}
}
\details{
Predictions are calculated from the \code{\link{grpnet}} object fit to the full sample of data, which is stored as \code{object$grpnet.fit} 

See \code{\link{predict.grpnet}} for further details on the calculation of the different types of predictions.
}
\value{
Depends on three factors... \cr
1. the exponential family distribution \cr
2. the length of the input \code{s} \cr
3. the \code{type} of prediction requested

See \code{\link{predict.grpnet}} for details
}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}

Helwig, N. E. (2024). Versatile descent algorithms for group regularization and variable selection in generalized linear models. \emph{Journal of Computational and Graphical Statistics}. \doi{10.1080/10618600.2024.2362232}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Syntax is inspired by the \code{predict.cv.glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010). 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}

\code{\link{predict.grpnet}} for predicting from \code{\link{grpnet}} objects
}
\examples{
\donttest{
######***######   family = "gaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = mpg)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto)

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto)

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(auto$mpg - fit.1se))
mean(abs(auto$mpg - fit.min))



######***######   family = "multigaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = (mpg, displacement))
y <- as.matrix(auto[,c(1,3)])
set.seed(1)
mod <- cv.grpnet(y ~ ., data = auto[,-c(1,3)], family = "multigaussian",
                 standardize.response = TRUE)

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto)

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(y - fit.1se))
mean(abs(y - fit.min))



######***######   family = "binomial"   ######***######

# load data
data(auto)

# redefine origin (Domestic vs Foreign)
auto$origin <- ifelse(auto$origin == "American", "Domestic", "Foreign")

# 10-fold cv (default method, response = origin with 2 levels)
set.seed(1)
mod <- cv.grpnet(origin ~ ., data = auto, family = "binomial")

# get predicted classes at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "class")

# get predicted classes at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "class", s = "lambda.min")

# compare misclassification rate for two solutions
1 - mean(auto$origin == fit.1se)
1 - mean(auto$origin == fit.min)



######***######   family = "multinomial"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin with 3 levels)
set.seed(1)
mod <- cv.grpnet(origin ~ ., data = auto, family = "multinomial")

# get predicted classes at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "class")

# get predicted classes at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "class", s = "lambda.min")

# compare misclassification rate for two solutions
1 - mean(auto$origin == fit.1se)
1 - mean(auto$origin == fit.min)



######***######   family = "poisson"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = horsepower)
set.seed(1)
mod <- cv.grpnet(horsepower ~ ., data = auto, family = "poisson")

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "response")

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "response", s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(auto$horsepower - fit.1se))
mean(abs(auto$horsepower - fit.min))



######***######   family = "negative.binomial"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = horsepower)
set.seed(1)
mod <- cv.grpnet(horsepower ~ ., data = auto, family = "negative.binomial")

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "response")

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "response", s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(auto$horsepower - fit.1se))
mean(abs(auto$horsepower - fit.min))



######***######   family = "Gamma"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto, family = "Gamma")

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "response")

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "response", s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(auto$mpg - fit.1se))
mean(abs(auto$mpg - fit.min))



######***######   family = "inverse.gaussian"   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = origin)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto, family = "inverse.gaussian")

# get fitted values at "lambda.1se"
fit.1se <- predict(mod, newdata = auto, type = "response")

# get fitted values at "lambda.min"
fit.min <- predict(mod, newdata = auto, type = "response", s = "lambda.min")

# compare mean absolute error for two solutions
mean(abs(auto$mpg - fit.1se))
mean(abs(auto$mpg - fit.min))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
