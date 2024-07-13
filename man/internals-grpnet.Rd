\name{internals-grpnet}

\alias{R_grpnet_binomial}
\alias{R_grpnet_gamma}
\alias{R_grpnet_gaussian}
\alias{R_grpnet_invgaus}
\alias{R_grpnet_multinom}
\alias{R_grpnet_negbin}
\alias{R_grpnet_poisson}

\alias{R_grpnet_binomial_dev}
\alias{R_grpnet_gamma_dev}
\alias{R_grpnet_invgaus_dev}
\alias{R_grpnet_multinom_dev}
\alias{R_grpnet_negbin_dev}
\alias{R_grpnet_poisson_dev}

\alias{R_grpnet_maxeigval}
\alias{R_grpnet_penalty}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal 'grpnet' Functions
}
\description{
Internal functions for 'grpnet' package, which are \code{R} translations of \code{Fortran} functions.
}
\usage{
# R translation of grpnet_binomial.f90
R_grpnet_binomial(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                  nlam, lambda, lmr, penid, gamma, eps, maxit,
                  standardize, intercept, ibeta, betas, iters,
                  nzgrps, nzcoef, edfs, devs, nulldev)

# R translation of grpnet_gamma.f90
R_grpnet_gamma(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
               nlam, lambda, lmr, penid, gamma, eps, maxit,
               standardize, intercept, ibeta, betas, iters,
               nzgrps, nzcoef, edfs, devs, nulldev)
             
# R translation of grpnet_gaussian.f90
R_grpnet_gaussian(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                  nlam, lambda, lmr, penid, gamma, eps, maxit,
                  standardize, intercept, ibeta, betas, iters,
                  nzgrps, nzcoef, edfs, devs, nulldev)             

# R translation of grpnet_invgaus.f90
R_grpnet_invgaus(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                 nlam, lambda, lmr, penid, gamma, eps, maxit,
                 standardize, intercept, ibeta, betas, iters,
                 nzgrps, nzcoef, edfs, devs, nulldev)

# R translation of grpnet_multinom.f90
R_grpnet_multinom(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw,  
                  alpha, nlam, lambda, lmr, penid, gamma, eps, maxit,
                  standardize, intercept, ibeta, betas, iters,
                  nzgrps, nzcoef, edfs, devs, nulldev)
                
# R translation of grpnet_negbin.f90
R_grpnet_negbin(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                nlam, lambda, lmr, penid, gamma, eps, maxit,
                standardize, intercept, ibeta, betas, iters,
                nzgrps, nzcoef, edfs, devs, nulldev, theta)                

# R translation of grpnet_poisson.f90
R_grpnet_poisson(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                 nlam, lambda, lmr, penid, gamma, eps, maxit,
                 standardize, intercept, ibeta, betas, iters,
                 nzgrps, nzcoef, edfs, devs, nulldev)                

# R translations of deviance functions
R_grpnet_binomial_dev(nobs, y, mu, wt, dev)
R_grpnet_gamma_dev(nobs, y, mu, wt, dev)
R_grpnet_invgaus_dev(nobs, y, mu, wt, dev)
R_grpnet_multinom_dev(nobs, nresp, y, mu, wt, dev)
R_grpnet_negbin_dev(nobs, y, mu, wt, theta, dev)
R_grpnet_poisson_dev(nobs, y, mu, wt, dev)

# R translations of utility functions
R_grpnet_maxeigval(A, N, MAXEV)
R_grpnet_penalty(znorm, penid, penone, pentwo, gamma, shrink)
}

\details{
  These functions are not intended to be called by the typical user.
}

\keyword{ internal }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
