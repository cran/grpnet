
!   grpnet_multigaus.f90 - group elastic net (multivariate gaussian)
!   Nathaniel E. Helwig (helwig@umn.edu)
!   Department of Psychology and School of Statistics
!   University of Minnesota
!   Date: 2025-01-17


! INPUTS/OUTPUTS
!   nobs = number of observations (N)
!   nvars = number of variables (P)
!   nresp = number of response dimensions for y (M)
!   x = grouped predictor matrix (N,P)
!       note: (x.1, x.2, ..., x.K) where x.k is (N,P.k) with sum(P.k) = P
!   y = response matrix (N,M)
!   w = observation weight vector (N)
!       note: if min(w) = max(w), then weights are ignored
!   off = offset matrix (N,M)
!   ngrps = number of groups (K)
!   gsize = number of coefs in each group (K)
!           note: (P.1, P.2, ..., P.K) with sum(P.k) = P
!   pw = penalty weight vector (K)   (outputs xsdev)
!        note: set pw(k) = 0 to leave k-th group unpenalized
!   alpha = weight for L1 and L2 penalities (alpha = 1 for lasso/mcp/scad, alpha = 0 for ridge)
!           note: setting alpha between (0,1) gives elastic net
!   nlam = number of lambdas
!   lambda = sequence of lambda values
!            note: if max(lambda) <=0, then lambda.max is data dependent
!   lmr = lambda minimum ratio (lambda.min = lmr * lambda.max)
!         note: unless lambda is provided, lambda.max is data dependent
!   penid = penalty id: 1 = lasso, 2 = mcp, 3 = scad
!   gamma = additional hyper-parameter for mcd and scad penalties
!           note: gamma > 1 for mcp and gamma > 2 for scad
!   eps = convergence tolerance
!   maxit = maximum number of iterations
!   standardize = integer (0: no standardization, 1: yes standardization)
!   intercept = integer (0: no intercept, 1: yes intercept)
!   ibeta = output matrix of intercepts (nresp,nlam)
!   betas = output array of coefficients (nvars,nresp,nlam)
!   iters = output vector of iterations for each lambda (nlam)
!   nzgrps = number of non-zero groups for each lambda (nlam)
!   nzcoef = number of non-zero coefficients for each lambda (nlam)
!   edfs = effective degrees of freedom for each lambda (nlam)
!   devs = residual deviance for each lambda (nlam)
!   nulldev = null deviance


SUBROUTINE grpnet_multigaus(nobs, nvars, nresp, x, y, w, off, ngrps, gsize, pw, &
                            alpha, nlam, lambda, lmr, penid, gamma, eps, maxit, &
                            standardize, intercept, ibeta, betas, iters, &
                            nzgrps, nzcoef, edfs, devs, nulldev)

    IMPLICIT NONE

! --------------- ARGUMENTS --------------- !
    INTEGER nobs, nvars, nresp, ngrps, gsize(ngrps), nlam, penid, maxit
    INTEGER standardize, intercept, iters(nlam), nzgrps(nlam), nzcoef(nlam)
    DOUBLE PRECISION x(nobs, nvars), y(nobs, nresp), w(nobs), off(nobs, nresp)
    DOUBLE PRECISION pw(ngrps), alpha, lambda(nlam), lmr, gamma, eps, ibeta(nresp,nlam)
    DOUBLE PRECISION betas(nvars, nresp, nlam), edfs(nlam), devs(nlam), nulldev
! --------------- ARGUMENTS --------------- !


! --------------- LOCAL DEFINITIONS --------------- !
    INTEGER i, j, k, l, iter, violations, makelambda
    INTEGER active(ngrps), strong(ngrps), ia(ngrps), ib(ngrps), gid
    DOUBLE PRECISION wmin, wmax, rnglam, minlam, maxlam, maxdif, macheps
    DOUBLE PRECISION r(nobs, nresp), difbeta(nvars, nresp), beta(nvars, nresp)
    DOUBLE PRECISION difibeta(nresp), zvec(nvars, nresp), grad(nvars, nresp)
    DOUBLE PRECISION gradnorm(ngrps), ctol, shrink, twolam, penone, pentwo
    DOUBLE PRECISION xmean(nvars), xsdev(ngrps), xev(ngrps), znorm, bnorm
    DOUBLE PRECISION wmat(nobs, nresp)
    DOUBLE PRECISION, ALLOCATABLE :: xtx(:,:)
! --------------- LOCAL DEFINITIONS --------------- !


! --------------- MACHINE EPSILON --------------- !
    macheps = EPSILON(eps)
! --------------- MACHINE EPSILON --------------- !


! --------------- CHECK WEIGHTS --------------- !
    wmin = MINVAL(w)
    wmax = MAXVAL(w)
    wmat = 1.0D0
    IF (wmax > wmin) THEN
        w = nobs * w / SUM(w)     ! normalize so SUM(w) = nobs
        w = SQRT(w)
        wmat = SPREAD(w, 2, nresp)
        y = wmat * y
        off = wmat * off
        DO i=1,nobs
            x(i,:) = w(i) * x(i,:)
        END DO
    END IF
! --------------- CHECK WEIGHTS --------------- !
    

! --------------- GROUP INDICES --------------- !
    gid = 0
    DO k=1,ngrps
        ia(k) = gid + 1
        gid = gid + gsize(k)
        ib(k) = gid
    END DO
! --------------- GROUP INDICES --------------- !


! --------------- CENTER AND SCALE --------------- !
    IF (intercept == 1) THEN
        DO j=1,nvars
            xmean(j) = SUM(x(:,j)) / nobs
            x(:,j) = x(:,j) - xmean(j)
        END DO
    END IF
    xsdev = 1.0D0
    IF (standardize == 1) THEN
        DO k=1,ngrps
            xsdev(k) = SQRT( SUM(x(:, ia(k):ib(k))**2) / (nobs * gsize(k)) )
            IF (xsdev(k) > macheps) THEN
                x(:,ia(k):ib(k)) = x(:,ia(k):ib(k)) / xsdev(k)
            ELSE
                xsdev(k) = 1.0D0
            END IF
        END DO
    END IF
! --------------- CENTER AND SCALE --------------- !


! --------------- GET MAX EIGENVALUE --------------- !
    DO k=1,ngrps
        IF (gsize(k) == 1) THEN
            xev(k) = SUM(x(:,ia(k))**2) / nobs
        ELSE
            ALLOCATE(xtx(gsize(k), gsize(k)))
            xtx = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), x(:,ia(k):ib(k))) / nobs
            CALL grpnet_maxeigval(xtx, gsize(k), xev(k))
            DEALLOCATE(xtx)
        END IF
    END DO
! --------------- GET MAX EIGENVALUE --------------- !


! --------------- MISCELLANEOUS INITIALIZATIONS --------------- !
    maxlam = MAXVAL(lambda)
    makelambda = 0
    iter = 0
    strong = 0
    active = 0
    nzgrps = 0
    nzcoef = 0
    zvec = 0.0D0
    ibeta = 0.0D0
    beta = 0.0D0
    grad = 0.0D0
    gradnorm = 0.0D0
    twolam = 0.0D0
    devs = 0.0D0
    r = y - off
! --------------- MISCELLANEOUS INITIALIZATIONS --------------- !


! --------------- GENERATE LAMBDA --------------- !
    IF (maxlam <= macheps) THEN

        makelambda = 1
        i = 1

        ! find unpenalized groups !
        DO k=1,ngrps
            IF (pw(k) <= macheps) THEN
                active(k) = 1
                nzgrps(i) = nzgrps(i) + 1
                nzcoef(i) = nzcoef(i) + gsize(k)
            END IF
        END DO
        ! find unpenalized groups !

        ! iterate until active coefficients converge !
        IF (nzgrps(i) > 0) THEN
            DO WHILE(iter < maxit)

                ! update iter and reset counters
                ctol = 0.0D0
                iter = iter + 1

                ! update active groups
                DO k=1,ngrps
                    IF(active(k) == 0) CYCLE
                    grad(ia(k):ib(k),:) = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), r) / nobs
                    zvec(ia(k):ib(k),:) = beta(ia(k):ib(k),:) + grad(ia(k):ib(k),:) / xev(k)
                    difbeta(ia(k):ib(k),:) = zvec(ia(k):ib(k),:) - beta(ia(k):ib(k),:)
                    maxdif = MAXVAL( ABS(difbeta(ia(k):ib(k),:)) / (1.0D0 + ABS(beta(ia(k):ib(k),:))) )
                    beta(ia(k):ib(k),:) = beta(ia(k):ib(k),:) + difbeta(ia(k):ib(k),:)
                    r = r - MATMUL(x(:,ia(k):ib(k)), difbeta(ia(k):ib(k),:))
                    ctol = MAX(maxdif , ctol)
                END DO ! k=1,ngrps

                ! update intercept
                IF (intercept == 1) THEN
                    difibeta = SUM(r * wmat, 1) / nobs
                    maxdif = MAXVAL( ABS(difibeta) / (1.0D0 + ABS(ibeta(:,i))) )
                    ibeta(:,i) = ibeta(:,i) + difibeta
                    r = r - wmat * SPREAD(difibeta, 1, nobs)
                    ctol = MAX(maxdif, ctol)
                END IF ! (intercept == 1)

                ! convergence check
                IF(ctol < eps) EXIT

            END DO ! WHILE(iter < maxit)

        ELSE

            ! intercept only
            iter = 1
            IF (intercept == 1) THEN
                ibeta(:,i) = SUM(r * wmat, 1) / nobs
                r = r - wmat * SPREAD(ibeta(:,i), 1, nobs)
            END IF

        END IF ! (nzgrps(1) > 0)
        ! iterate until active coefficients converge !

        ! create lambda sequence !
        DO k=1,ngrps
            IF (pw(k) > macheps) THEN
                grad(ia(k):ib(k),:) = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), r) / nobs
                gradnorm(k) = SQRT(SUM(grad(ia(k):ib(k),:)**2)) / pw(k)
            END IF
        END DO
        IF (alpha > macheps) THEN
            maxlam = MAXVAL(gradnorm / alpha)
        ELSE
            maxlam = MAXVAL(gradnorm / 1.0E-3)
            makelambda = 0
        END IF
        minlam = lmr * maxlam
        lambda(i) = maxlam
        maxlam = LOG(maxlam)
        minlam = LOG(minlam)
        rnglam = maxlam - minlam
        DO k=2,nlam
            lambda(k) = EXP(maxlam - rnglam * (k - 1) / (nlam - 1))
        END DO
        ! create lambda sequence !

        ! save results !
        betas(:,:,i) = beta
        iters(i) = iter
        nzgrps(i) = nzgrps(i) + intercept
        nzcoef(i) = nzcoef(i) + intercept * nresp
        edfs(i) = DBLE(nzcoef(i))
        devs(i) = SUM(r**2)
        ! save results !

    END IF
! --------------- GENERATE LAMBDA --------------- !


! --------------- ITERATIVE WORK --------------- !
    DO i=1,nlam

        ! initializations !
        IF (i == 1 .AND. makelambda == 1) CYCLE
        IF (i > 1) THEN
            ibeta(:,i) = ibeta(:,i-1)
            beta = betas(:,:,i-1)
            twolam = alpha * (2.0D0 * lambda(i) - lambda(i-1))
        ELSE
            grad = MATMUL(TRANSPOSE(x), r) / nobs
        END IF
        ! initializations !

        ! strong rule initialization !
        DO k=1,ngrps
            gradnorm(k) = SQRT(SUM(grad(ia(k):ib(k),:)**2))
            IF (gradnorm(k) + 1.0E-8 > pw(k) * twolam) THEN
                strong(k) = 1
            ELSE
                strong(k) = 0
            END IF
        END DO
        ! strong rule initialization !

        ! iterate until strong set converges !
        iter = 0
        DO WHILE(iter < maxit)

            ! iterate until active set converges !
            DO WHILE(iter < maxit)

                ! iterate until active coefficients converge !
                DO WHILE(iter < maxit)

                    ! update iter and reset counters
                    iter = iter + 1
                    nzgrps(i) = 0
                    edfs(i) = 0.0D0
                    ctol = 0.0D0

                    ! update active groups
                    DO k=1,ngrps
                        IF(active(k) == 0) CYCLE
                        penone = alpha * lambda(i) * pw(k) / xev(k)
                        pentwo = (1.0D0 - alpha) * lambda(i) * pw(k) / xev(k)
                        grad(ia(k):ib(k),:) = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), r) / nobs
                        zvec(ia(k):ib(k),:) = beta(ia(k):ib(k),:) + grad(ia(k):ib(k),:) / xev(k)
                        znorm = SQRT(SUM(zvec(ia(k):ib(k),:)**2))
                        bnorm = SQRT(SUM(beta(ia(k):ib(k),:)**2))
                        CALL grpnet_penalty(znorm, penid, penone, pentwo, gamma, shrink)
                        IF(shrink == 0.0D0 .AND. bnorm == 0.0D0) CYCLE
                        difbeta(ia(k):ib(k),:) = shrink * zvec(ia(k):ib(k),:) - beta(ia(k):ib(k),:)
                        maxdif = MAXVAL( ABS(difbeta(ia(k):ib(k),:)) / (1.0D0 + ABS(beta(ia(k):ib(k),:))) )
                        beta(ia(k):ib(k),:) = beta(ia(k):ib(k),:) + difbeta(ia(k):ib(k),:)
                        r = r - MATMUL(x(:,ia(k):ib(k)), difbeta(ia(k):ib(k),:))
                        ctol = MAX(maxdif , ctol)
                        IF(shrink > 0.0D0) THEN
                            nzgrps(i) = nzgrps(i) + 1
                            edfs(i) = edfs(i) + gsize(k) * shrink
                        ENDIF
                    END DO ! k=1,ngrps

                    ! update intercept
                    IF (intercept == 1) THEN
                        difibeta = SUM(r * wmat, 1) / nobs
                        maxdif = MAXVAL( ABS(difibeta) / (1.0D0 + ABS(ibeta(:,i))) )
                        ibeta(:,i) = ibeta(:,i) + difibeta
                        r = r - wmat * SPREAD(difibeta, 1, nobs)
                        ctol = MAX(maxdif, ctol)
                    END IF ! (intercept == 1)

                    ! convergence check
                    IF(ctol < eps) EXIT

                END DO ! WHILE(iter < maxit) - inner
                ! iterate until active coefficients converge !

                ! check inactive groups in strong set !
                violations = 0
                DO k=1,ngrps
                    IF(strong(k) == 0 .OR. active(k) == 1) CYCLE
                    grad(ia(k):ib(k),:) = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), r) / nobs
                    gradnorm(k) = SQRT(SUM(grad(ia(k):ib(k),:)**2))
                    IF (gradnorm(k) > alpha * lambda(i) * pw(k)) THEN
                        active(k) = 1
                        violations = violations + 1
                    END IF
                END DO ! k=1,ngrps
                IF(violations == 0) EXIT
                ! check inactive groups in strong set !

            END DO ! WHILE(iter < maxit) - middle
            ! iterate until active set converges !

            ! check groups in weak set !
            violations = 0
            DO k=1,ngrps
                IF (strong(k) == 1) CYCLE
                grad(ia(k):ib(k),:) = MATMUL(TRANSPOSE(x(:,ia(k):ib(k))), r) / nobs
                gradnorm(k) = SQRT(SUM(grad(ia(k):ib(k),:)**2))
                IF (gradnorm(k) + 1e-8 > alpha * lambda(i) * pw(k)) THEN
                    strong(k) = 1
                    active(k) = 1
                    violations = violations + 1
                END IF
            END DO
            IF(violations == 0) EXIT
            ! check groups in weak set !

        END DO ! WHILE(iter < maxit) - outer
        ! iterate until strong set converges !

        ! calculate nzcoef !
        DO k=1,ngrps
            IF(active(k) == 0) CYCLE
            DO j=1,nresp
                DO l=1,gsize(k)
                    IF (ABS(beta(ia(k) + l - 1,j)) > macheps) THEN
                        nzcoef(i) = nzcoef(i) + 1
                    END IF
                END DO
            END DO
        END DO
        ! calculate nzcoef !

        ! save results !
        betas(:,:,i) = beta
        iters(i) = iter
        nzgrps(i) = nzgrps(i) + intercept
        nzcoef(i) = nzcoef(i) + intercept * nresp
        edfs(i) = edfs(i) + DBLE(intercept * nresp)
        devs(i) = SUM(r**2)
        ! save results !

    END DO
! --------------- ITERATIVE WORK --------------- !


! --------------- POST PROCESSING --------------- !
    IF (standardize == 1) THEN
        DO k=1,ngrps
            betas(ia(k):ib(k),:,:) = betas(ia(k):ib(k),:,:) / xsdev(k)
        END DO
    END IF
    IF (intercept == 1) THEN
        DO j=1,nresp
            ibeta(j,:) = ibeta(j,:) - MATMUL(xmean, betas(:,j,:))
        END DO
        off = wmat * SPREAD(SUM(y * wmat, 1) / nobs, 1, nobs)
        nulldev = SUM((y - off)**2)
    ELSE
        nulldev = SUM((y - off)**2)
    END IF
    pw = xsdev
! --------------- POST PROCESSING --------------- !


END SUBROUTINE
