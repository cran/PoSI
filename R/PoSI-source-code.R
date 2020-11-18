#================================================================
#
#
#      R SOFTWARE FOR COMPUTING POSI CONSTANTS: SOURCE CODE
#
#
# Authors: Andreas Buja, Kai Zhang
# 2014/12/11
#
# Reference: http://arxiv.org/abs/1306.1059
#            "Valid Post-Selection Inference,"
#            by Berk, R., Brown, L., Buja, A., Zhang, K., Zhao,L.,
#            The Annals of Statistics, 41 (2), 802-837 (2013)
#
#
#================================================================
##
##
## PoSI
##
## Valid Post-Selection Inference for multiple linear regression
##
## Description:
##
##     'PoSI()' creates an object of class "PoSI" from a design/predictor
##     matrix.  PoSI objects are used to obtain the standard error
##     multipliers or z-test/t-test thresholds for valid post-selection
##     inference in multiple linear regression according to Berk et
##     al. (2013, Annals of Statistics).  These multipliers are valid
##     after performing arbitrary variable selection in the specified set
##     of submodels.
##
##     The functions 'summary()' and 'K()' compute the actual multipliers
##     from an object of class 'PoSI' returned by the function 'PoSI()'.
##     The function 'summary()' provides PoSI multipliers for given
##     confidence levels as well as more conservative multipliers based
##     on the Scheffe method of simultaneous inference and Bonferroni
##     adjustment; 'K()' provides just the PoSI multipliers.  By default
##     the multipliers are provided for z-tests assuming the error
##     variance is known; if error degrees of freedom are specified
##     ('df.err'), the multipliers are for t-tests.
##
##
## Usage:
##
##     PoSI(X, modelSZ=1:ncol(X), Nsim=1000, bundleSZ=10000,
##          eps=1E-8, center=T, scale=T, verbose=1)
##
##     summary(PoSI.object, confidence=c(.95,.99), alpha=NULL, df.err=NULL,
##             eps.PoSI=1E-6, digits=3)
##
##     K(PoSI.object, confidence=c(.95,.99), alpha=NULL, df.err=NULL,
##       eps.PoSI=1E-6, digits=3
##
## Arguments:
##
##            X: a design/predictor/covariate matrix for multiple linear
##               regression
##
##               [When there are predictors that are 'protected', that
##               is, forced into the model and not subjected to variable
##               selection, one should adjust the selectable predictors
##               for the protected predictors (i.e., replace the
##               selectable predictors with their residual vectors when
##               regressing them on the protected predictors); thereafter
##               simply omit the protected predictors from X.]
##
##      modelSZ: vector of submodel sizes to be searched
##               Default: 1:ncol(X), i.e., all submodels of all sizes.
##               If only models up to size 5 (e.g.) are considered, use 1:5.
##
##         Nsim: the number of random unit vectors sampled in the column
##               space of X; default: Nsim=1000
##               'PoSI' is simulation-based. Larger 'Nsim' increases
##               precision but may cause running out of memory.
##
##     bundleSZ: the number of contrasts/adjusted predictor vectors
##               processed at any given time; default: bundleSZ=10000
##               'PoSI' is more efficient with larger 'bundleSZ', but it
##               may cause running out of memory.
##
##          eps: threshold below which singular values of X are
##               considered to be zero
##               In case of highly collinear columns of X this threshold
##               determines the effective dimension of the column space
##               of X.  Default: eps=1E-6
##
##       center: whether the columns of X should be centered
##               In most applications the intercept is not among
##               variables to be selected.  Centering effectively adjusts
##               for the presence of an intercept that should be included
##               in all submodels.
##
##        scale: whether the columns of X should be standardized
##               This is appropriate to prevent numeric problems from
##               predictors with vastly differing scales.
##
##      verbose: verbosity reporting levels
##               0: no reports
##               1: report bundle completion (default);
##               2: report each processed submodel (for debugging)
##
##   confidence: list of confidence levels for which standard error
##               multipliers should be provided
##               default: confidence=c(.95,.99)
##
##        alpha: if specified, sets  'confidence = 1-alpha'
##
##       df.err: error degrees of freedom for t-tests
##               default: df.err=NULL (z-tests)
##
##     eps.PoSI: precision to which standard error multipliers are
##               computed
##
##       digits: number of digits to which standard error multipliers are
##               rounded
##
## Details:
##
##       'PoSI' is simulation-based and requires specification of a
##       number 'Nsim' of random unit vectors to be sampled in the column
##       space of 'X'.  Large 'Nsim' yields greater precision but
##       requires more memory.  The memory demands can be lowered by
##       decreasing 'bundleSZ' at the cost of some efficiency.
##       'bundleSZ' determines how many adjusted predictor vectors are
##       being processed at any given time.
##
##       The computational limit of the PoSI method is in the exponential
##       growth of the number of adjusted predictor vectors:
##       - If p=ncol(X) and all submodels are being searched, the number
##         of adjusted predictor vectors is p*2^(p-1).
##         Example: p=20, m=1:20  ==>  |{adj.vecs.}| = 10,485,760
##       - If only models of a given size 'm' are being searched, the number
##         of adjusted predictor vectors is m*choose(p,m).
##         Example: p=50, m=1:5   ==>  |{adj.vecs.}| = 11,576,300
##       Thus limiting PoSI to small submodel sizes such as 'm=1:4'
##       ("sparse PoSI") puts larger p=ncol(X) within reach.
##
## Value:
##
##       'PoSI' returns an object of 'class' '"PoSI"'.
##
##       'summary' returns a small matrix with standard error multipliers
##       for all requested confidence levels as well as computed by three
##       different methods: the PoSI method, Scheffe simultaneous
##       inference, and Bonferroni adjustment.
##
##       'K' returns just the PoSI-based standard error multipliers for
##       the requested confidence levels or for specified Type I errors
##       'alpha' (default: NULL, overwrites 'confidence = 1 - alpha' if
##       non-NULL).
##
##
#================================================================


# Computations of inner products of random unit vectors (U)
# with 'linear contrasts', that is, adjusted predictor vectors (L):
PoSI <- function(X, modelSZ=1:ncol(X), center=T, scale=T, verbose=1,
                 Nsim=1000, bundleSZ=100000, eps=1E-8)
{
    if(verbose>=1) { proc.time.init <- proc.time() }
    # Standardize the columns of X for numerical stability; if intercept is to be selected, set 'center=F':
    X.scaled <- scale(X, center=center, scale=scale)
    # Remove columns with NAs:
    sel <- apply(X.scaled, 2, function(x) !any(is.na(x)))
    X.scaled <- X.scaled[,sel]
    if(verbose >= 1 & any(!sel)) cat("Removed column(s)",which(!sel),"due to NAs after standardizing'\n")
    # Map X to Scheffe canonical coordinates X.cc:  X.cc' X.cc = X'X  (d x p)
    X.svd   <- svd(X.scaled)
    sel <- (X.svd$d > eps)
    if(verbose >= 1 & any(!sel)) cat("Reduced rank:",sum(sel),"  (Change eps?)\n")
    X.cc <- X.svd$d[sel] * t(X.svd$v[,sel]) # Full-model Scheffe canonical coords. of X.
    # Clean up:
    rm(X.scaled, X.svd)
    # Notation:
    p <- ncol(X.cc)                         # Number of predictors
    d <- nrow(X.cc)                         # d=rank(X), can be < p when X is rank-deficient, as when p>n
    # Submodels can only be of size <= d; drop those exceeding d:
    modelSZ <- modelSZ[modelSZ <= d]        # Remove models that are too large.
    # Random unit vectors for which inner products with adjusted predictors (L) will be formed:
    U <- matrix(rnorm(d*Nsim),ncol=d)       # Nsim x d
    U <- U / sqrt(rowSums(U^2))             # Matrix of unit vectors u in R^d (=rows) for max_l |<u,l>|
    # Prepare loop over models to store all possible linear combinations l:
    UL.max <- rep(0,Nsim)                   # To contain  max_l |<u,l>|  for each random row u of U.
    L <- matrix(NA, nrow=bundleSZ, ncol=d)  # bundleSZ x d, one bundle of adjusted predictors l
    Nmodels  <- 0                           # Number of processed models; to be returned among attributes()
    Ntests   <- 0                           # Number of processed t-tests; to be returned among attributes()
    ibundles <- 0                           # Bundle counter, to be reported when 'verbose >= 1'.
    Nbundles <- sum(ceiling(choose(p,modelSZ)*modelSZ/bundleSZ))
    if(verbose>=1) {
        cat("Number of contrasts/adjusted predictors to process:",sum(choose(p,modelSZ)*modelSZ),"\n")
        cat("Number of bundles:",Nbundles,"\n")
    }
    for(m in modelSZ) {                     # Loop over model sizes.
        if(m > bundleSZ) {                  # Silly: bundleSZ is too small to contain a single model; fix it.
            bundleSZ <- max(bundleSZ, m)
            L <- matrix(NA, nrow=bundleSZ, ncol=d)
        }
        i.store <- 1                        # Pointer into next free row of L (= d-vectors)
        M <- 1:m                            # M = current submodel, initialized to first submodel
        M.last <- (p-m+1):p                 # M.last = last submodel of size m
        repeat {                            # Loop over bundles of adjusted predictors
            # Verbosity level 2: print all models & dims
            if(verbose>=2) cat("Model:", M,"\n")
            X.m <- X.cc[,M]                 # Select predictors for this submodel
            X.qr <- qr(t(X.m)%*%X.m)        # QR decomposition of submodel in canon. coords.
            if(X.qr$rank == m) {            # Process only if full rank; else skip to next.
                L[i.store:(i.store+m-1),] <- t(X.m %*% qr.solve(X.qr)) # Submodel adjusted predictors (rows of L)
                Nmodels <- Nmodels + 1
                Ntests  <- Ntests + m
                i.store <- i.store+m        # Move pointer to next free row of L.
                # If bundle L of adjusted predictors l is full, or if this was the last model, then
                # accumulate  max_l |<u,l>|  and re-initialize L:
                if(( (i.store-1+m) > bundleSZ ) | !(any(M < M.last) )) {
                    ii.store <- 1:(i.store-1)                     # Filled rows in L
                    L[ii.store,] <- L[ii.store,,drop=F] /         # Normalize rows of L to unit length:
                        sqrt(rowSums(L[ii.store,,drop=F]^2))
                    # Implementation note: Against 'help()', '..%*%t(..)' is a touch faster than 'tcrossprod()'.
                    UL <- abs(U %*% t(L[ii.store,,drop=F]))         # Inner prods |<u,l>| of rows of U and L
                    # UL <- abs(tcrossprod(U, L[ii.store,,drop=F])) # Inner prods |<u,l>| of rows of U and L
                    UL.max <- pmax(UL.max, apply(UL, 1, max))     # Accumulate  max_l |<u,l>|
                    i.store <- 1                                  # Initialize L to fill up a new bundle
                    # Verbosity level 1: print on bundle completion
                    if(verbose>=1) {
                        ibundles <- ibundles + 1
                        cat("                         Done with bundle",ibundles,"/",Nbundles,"   model sz =",m,"\n")
                    }
                } # end of 'if(( (i.store-1+m...'
            } # end of 'if(X.qr$rank...'
            # Find next model in enumeration:
            if(any(M < M.last)) {            # Not at the last submodel yet, hence find the next one:
                i <- max(which(M < M.last))  # Pointer to the highest predictor index not in the last model
                M[i] <- M[i]+1   # Move this predictor index to point to the next predictor.
                if(i < m) { M[(i+1):m] <- (M[i]+1):(M[i]+m-i) } # Set the remaining to their minima.
            } else break
        } # end of 'repeat{'
    } # end of 'for(m in modelSZ) {'
    if(verbose>=1) {
        cat("p =",p,", d =",d,"  processed",Ntests,"tests in",Nmodels,"models.  Times in seconds:\n")
        print(proc.time() - proc.time.init)
    }
    attributes(UL.max)$d <- d
    attributes(UL.max)$p <- p
    attributes(UL.max)$Nmodels <- Nmodels
    attributes(UL.max)$Ntests  <- Ntests
    class(UL.max) <- "PoSI"
    UL.max
} # end of 'PoSI <- function(...)'


#================================================================


# Computation of PoSI, Scheffe and Bonferroni constants:
summary.PoSI <- function(object, confidence=c(.95,.99), alpha=NULL, df.err=NULL,
                         eps.PoSI=1E-6, digits=3, ...)
{
    if(class(object) != "PoSI") {
        cat("summary.PoSI: first argument is not of class 'PoSI'.\n")
        return()
    }
    rU <- 1/object
    df <- attributes(object)$d
    if(is.null(alpha)) { alpha <- 1-confidence } else { confidence <- 1-alpha }
    if(any(alpha >= 1) | any(confidence >= 1)) {
        cat("summary.PoSI: args 'confidence' or 'alpha' out of range (0,1).\n")
        return()
    }
    # The following functions rely on lexical scoping:
    if(is.null(df.err)) { # sigma known  ==>  z-statistics
        prob.PoSI  <- function(K) mean(pchisq((K*rU)^2, df))
        quant.Sch  <- function() sqrt(qchisq(confidence, df=df))
        quant.Bonf <- function() qnorm(1 - alpha/2/attributes(object)$Ntests)
    } else { # sigma estimated with df.err  ==>  t-statistics
        prob.PoSI  <- function(K) mean(pf((K*rU)^2/df, df, df.err))
        quant.Sch  <- function() sqrt(df*qf(confidence, df1=df, df2=df.err))
        quant.Bonf <- function() qt(1 - alpha/2/attributes(object)$Ntests, df=df.err)
    }
    # Bisection search for K.PoSI:
    K.PoSI <- sapply(confidence,
                     function(cvg) {
                         K.lo <- 0
                         K.hi <- 30
                         while(abs(K.lo - K.hi) > eps.PoSI) {
                             K <- (K.lo + K.hi)/2
                             if(mean(prob.PoSI(K)) > cvg) {
                                 K.hi <- K
                             } else {
                                 K.lo <- K
                             } }
                         K }
                     )
    K.Sch <- quant.Sch()
    K.Bonf <- quant.Bonf()
    Ks <- cbind("K.PoSI"=K.PoSI, "K.Bonferroni"=K.Bonf, "K.Scheffe"=K.Sch)
    rownames(Ks) <- paste(confidence*100,"%",sep="")
    round(Ks, digits)
} # end of 'summary.PoSI <- function(...'

#================================================================
