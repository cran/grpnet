\name{internals-grpnet}

\alias{grpnet_binomial}
\alias{grpnet_binomial_dev}
\alias{grpnet_gamma}
\alias{grpnet_gamma_dev}
\alias{grpnet_gaussian}
\alias{grpnet_invgaus}
\alias{grpnet_invgaus_dev}
\alias{grpnet_logit}
\alias{grpnet_logit_dev}
\alias{grpnet_maxeigval}
\alias{grpnet_multigaus}
\alias{grpnet_multinom}
\alias{grpnet_multinom_dev}
\alias{grpnet_negbin}
\alias{grpnet_negbin_dev}
\alias{grpnet_ordinal}
\alias{grpnet_ordinal_int}
\alias{grpnet_penalty}
\alias{grpnet_poisson}
\alias{grpnet_poisson_dev}
\alias{grpnet_svm1}
\alias{grpnet_svm1_dev}
\alias{grpnet_svm2}
\alias{grpnet_svm2_dev}

\alias{R_grpnet_binomial}
\alias{R_grpnet_binomial_dev}
\alias{R_grpnet_gamma}
\alias{R_grpnet_gamma_dev}
\alias{R_grpnet_gaussian}
\alias{R_grpnet_invgaus}
\alias{R_grpnet_invgaus_dev}
\alias{R_grpnet_logit}
\alias{R_grpnet_logit_dev}
\alias{R_grpnet_maxeigval}
\alias{R_grpnet_multigaus}
\alias{R_grpnet_multinom}
\alias{R_grpnet_multinom_dev}
\alias{R_grpnet_negbin}
\alias{R_grpnet_negbin_dev}
\alias{R_grpnet_penalty}
\alias{R_grpnet_ordinal}
\alias{R_grpnet_ordinal_int}
\alias{R_grpnet_poisson}
\alias{R_grpnet_poisson_dev}
\alias{R_grpnet_svm1}
\alias{R_grpnet_svm1_dev}
\alias{R_grpnet_svm2}
\alias{R_grpnet_svm2_dev}


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
                 
# R translation of grpnet_logit.f90
R_grpnet_logit(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
               nlam, lambda, lmr, penid, gamma, eps, maxit,
               standardize, intercept, ibeta, betas, iters,
               nzgrps, nzcoef, edfs, devs, nulldev)

# R translation of grpnet_multigaus.f90
R_grpnet_multigaus(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw, 
                   alpha, nlam, lambda, lmr, penid, gamma, eps, maxit,
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

# R translation of grpnet_ordinal.f90
R_grpnet_ordinal(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw,
                 alpha, nlam, lambda, lmr, penid, gamma, eps, maxit,
                 standardize, intercept, ibeta, betas, iters,
                 nzgrps, nzcoef, edfs, devs, nulldev)
                 
# R translation of grpnet_ordinal_int.f90
R_grpnet_ordinal_int(nobs, nresp, eta, y, w, yfreq, theta, maxdif)                 

# R translation of grpnet_poisson.f90
R_grpnet_poisson(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
                 nlam, lambda, lmr, penid, gamma, eps, maxit,
                 standardize, intercept, ibeta, betas, iters,
                 nzgrps, nzcoef, edfs, devs, nulldev)
                 
# R translation of grpnet_svm1.f90
R_grpnet_svm1(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
              nlam, lambda, lmr, penid, gamma, eps, maxit,
              standardize, intercept, ibeta, betas, iters,
              nzgrps, nzcoef, edfs, devs, nulldev, theta)
              
# R translation of grpnet_svm2.f90
R_grpnet_svm2(nobs, nvars, x, y, w, off, ngrps, gsize, pw, alpha, 
              nlam, lambda, lmr, penid, gamma, eps, maxit,
              standardize, intercept, ibeta, betas, iters,
              nzgrps, nzcoef, edfs, devs, nulldev)

# R translations of deviance functions
R_grpnet_binomial_dev(nobs, y, mu, wt, dev)
R_grpnet_gamma_dev(nobs, y, mu, wt, dev)
R_grpnet_invgaus_dev(nobs, y, mu, wt, dev)
R_grpnet_logit_dev(nobs, y, mu, wt, dev)
R_grpnet_multinom_dev(nobs, nresp, y, mu, wt, dev)
R_grpnet_negbin_dev(nobs, y, mu, wt, theta, dev)
R_grpnet_poisson_dev(nobs, y, mu, wt, dev)
R_grpnet_svm1_dev(nobs, y, mu, wt, theta, dev)
R_grpnet_svm2_dev(nobs, y, mu, wt, dev)

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
