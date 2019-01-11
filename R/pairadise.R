#' pairadise
#'
#'
#' Primary function of the PAIRADISE package. Analyzes matched pairs for differences in isoform expression.
#' Uses parallel processing to speed up computation.
#'
#' @name pairadise
#' @param pdat A PDseDataSet object
#' @param nIter Positive integer. Specifies the maximum number of iterations of the optimization
#' algorithm allowed. Default is nIter = 100
#' @param tol Positive number. Specifies the tolerance level for terminating the optimization
#' algorithm, defined as the difference in log-likelihood ratios between iterations. Default
#' is tol = 10^(-2)
#' @param pseudocount Positive number. Specifies a value for a pseudocount added to each
#' count at the beginning of the analysis. Default is pseudocount = 0
#' @param seed An integer to set seed.
#' @param equal.variance Are the group variances assumed equal? Default value is FALSE.
#' @param numCluster Number of clusters to use for parallel computing.
#' @param BPPARAM parallel parameters from package BiocParallel.
#' @details This is the primary function of the PAIRADISE package that implements the PAIRADISE algorithm.
#' @return A PDseDataSet object contains outputs from PAIRADISE algorithm.
#' @examples
#'
#' #############################
#' ## Example: Simulated data ##
#' #############################
#'
#' set.seed(12345)
#' data("sample_dataset")
#' pdat <- PDseDataSetFromMat(sample_dataset)
#' pdat <- pairadise(pdat, numCluster =4)
#' results(pdat)
#' @export
#' @import S4Vectors
#' @import nloptr
#' @import BiocParallel
#' @importFrom methods is
#' @importFrom stats optim p.adjust pchisq rnorm
#' @rdname pairadise-FUN
pairadise <- function(pdat, nIter = 100, tol = 10^(-2), pseudocount = 0,
                      seed = 12321, equal.variance = FALSE, numCluster = 2,
                      BPPARAM = MulticoreParam(numCluster)) {
    stopifnot(is(pdat, "PDseDataSet"))
    
    if (nIter <= 0) {
        stop("Error: Number of iterations must be at least 1")
    }
    
    if (tol <= 0) {
        stop("Error: Tolerance must be strictly positive")
    }
    
    if (pseudocount < 0) {
        stop("Error: Pseudocount must be nonnegative")
    }
    
    if (!(is.logical(equal.variance))) {
        stop("Error: equal.variance must either be TRUE or FALSE")
    }

    if(numCluster > 1){
        outs <- bplapply(
            pdat, .pairadise, BPPARAM = BPPARAM,
            nIter = nIter, tol = tol,
            pseudocount = pseudocount, seed = seed,
            equal.variance = equal.variance)
    }else{
        outs <- lapply(
            pdat, .pairadise, nIter = nIter, tol = tol,
            pseudocount = pseudocount, seed = seed,
            equal.variance = equal.variance)
    }
    outs <- do.call(rbind, outs)
    rowData(pdat)$outs <- outs
    return(pdat)
}

#' Extract results for pairadise analysis
#' @param pdat A PDseDataSet object from pairadise analysis
#' @param p.adj The p ajustment method.
#' @param sig.level The cutoff of significant results
#' @param details Whether to list detailed results.
#' @return The function return a results DataFrame.
#'     \item{testStats}{Vector of test statistics for paired
#'     analysis.}  \item{p.value}{Vector of pvalues for each
#'     exon/event.}  \item{p.adj}{The adjusted p values} If details is
#'     TRUE, more detailed parameter estimates for constrained and
#'     unconstrained model will return.
#' @export
#' @examples
#' data("sample_dataset")
#' pdat <- PDseDataSetFromMat(sample_dataset)
#' pdat <- pairadise(pdat)
#' results(pdat)
results <- function(pdat, p.adj = "BH", sig.level = 0.01, details = FALSE){
    stopifnot(is(pdat, "PDseDataSet"))
    if (sig.level < 0 | sig.level > 1) {
        stop("Error: Significance level must be strictly between 0 and 1")
    }
    if(!"outs" %in% colnames(rowData(pdat))) stop("Please run pairadise first")
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~#~ Find significant exons ~#~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##    
    outs <- rowData(pdat)$outs

    paired.testStats <- unlist(outs[,1])
    ## All of the p-values of paired test
    total.pvals.paired <- 1 - pchisq(paired.testStats, 1)
    stat <- DataFrame(
        testStats=paired.testStats,
        p.value=total.pvals.paired)

    if(details){
        params <- data.frame(outs[,c(3,5:7,2, 4,8:10, 18)])
        colnames(params) <- c(
            "mu.u", "s1.u", "s2.u", "s.u", "delta",
            "mu.c", "s1.c", "s2.c", "s.c", "totalIter")

        latent <- lapply(seq_len(nrow(outs)), function(i){
            x <- outs[i,]
            lat <- do.call(rbind, (x[c(11:16)]))
            rownames(lat) <- paste(
                c("psi1", "psi2", "alpha"),
                rep(c("u", "c"), each=3), sep=".")
            lat
        })
        params$latent <- cbind(latent)
        stat <- DataFrame(stat, params)
    }

    if(!is.null(sig.level)){
        if(!is.null(p.adj)){
            stat$p.adj <- p.adjust(stat$p.value, method = p.adj)
            stat <- stat[stat$p.adj < sig.level,]
        }else{
            stat <- stat[stat$p.value < sig.level,]
        }
    }
    return(stat)
}

.pairadise <- function(pdat1, seed = 12321, nIter = 100, tol = 10^(-2),
                       pseudocount = 0, equal.variance = FALSE){

    l.iI <- rowData(pdat1)$iLen
    l.iS <- rowData(pdat1)$sLen
    count1 <- counts(pdat1)
    idx <- tapply(seq(ncol(pdat1)), colData(pdat1)$sample, c)
    idx1 <- do.call(rbind, idx)[,1]
    idx2 <- do.call(rbind, idx)[,2]
    count1 <- data.frame(cbind(count1[,idx1,], count1[,idx2,]))
    colnames(count1) <- c("I1", "S1", "I2", "S2")
    ## valid
    count1 <- na.omit(count1)
    count1 <- count1[(rowSums(count1[,1:2]) > 0) & (rowSums(count1[,3:4]) > 0),]
    
    I1 <- count1$I1 + pseudocount
    S1 <- count1$S1 + pseudocount
    I2 <- count1$I2 + pseudocount
    S2 <- count1$S2 + pseudocount
    M <- nrow(count1)
            
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~#~ Step 1: Initialize parameters ~#~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            
    ## Set seed
    set.seed(seed)
            
    ## The suffix .c always stands for 'constrained' The suffix .u always stands for
    ## 'unconstrained'
    logit.psi1.u <- logit.psi2.u <- logit.psi1.c <- logit.psi2.c <- 
        matrix(0, nrow = (nIter + 1), ncol = M)
    alpha.u <- alpha.c <- matrix(0, nrow = (nIter + 1), ncol = M)
            
    ## eps to prevent logit(0) or logit(1)
    eps <- 0.05
    logit.psi1.u[1, ] <- logit((I1 * l.iS + eps)/(I1 * l.iS + S1 * l.iI + 
                                                  2 * eps))
    logit.psi2.u[1, ] <- logit((I2 * l.iS + eps)/(I2 * l.iS + S2 * l.iI + 
                                                  2 * eps))
    alpha.u[1, ] <- logit.psi1.u[1, ] + rnorm(M, mean = 0, sd = 0.01)
            
    logit.psi1.c.old <- logit.psi1.u.old <- logit.psi1.c[1, ] <- logit.psi1.u[1, ]
    logit.psi2.c.old <- logit.psi2.u.old <- logit.psi2.c[1, ] <- logit.psi2.u[1, ]
    alpha.c.old <- alpha.u.old <- alpha.c[1, ] <- alpha.u[1, ]
            
    s1.u <- s1.c <- s2.u <- s2.c <- s.u <- s.c <- 0.1
    mu.u <- mu.c <- delta.u <- delta.c <- 0

    s1.u.old <- s1.c.old <- s.c.old <- s2.u.old <- s2.c.old <- s.u.old <- 0.1
    mu.u.old <- mu.c.old <- delta.u.old <- delta.c.old <- 0
            
    ll.old.u <- 10^(40)
    ll.old.c <- 10^(40)
            
    ## Define limits of parameters
    s.lower <- 10^(-6)
    s.upper <- Inf
    mu.lower <- -Inf
    mu.upper <- Inf
            
    MLE1.lower.u <- c(s.lower, s.lower, s.lower, mu.lower, -Inf)
    MLE1.lower.c <- c(s.lower, s.lower, s.lower, mu.lower, 0)
            
    MLE1.upper.u <- c(s.upper, s.upper, s.upper, mu.upper, Inf)
    MLE1.upper.c <- c(s.upper, s.upper, s.upper, mu.upper, 0)
            
    for (t in seq_len(nIter)) {
                
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~#~ Step 2: Estimate delta, mu, and sigmas ~#~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        
        ## This is the second step in the optimization process where we estimate the MLEs
        ## of delta, mu, sigma, sigma1, sigma2 based on the MLEs of logit(psi1),
        ## logit(psi2), and alpha computed in the previous stage.  See optimize1 for more
        ## details
        
        ## .u's always referred to unconstrained MLE .c's always referred to constrained
        ## MLE
                
        MLE1.u <- function(x) {
            optimize1(x, M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.u[t, ], 
                      logit.psi2.u[t, ], alpha.u[t, ], equal.variance)
        }
                
        MLE1.c <- function(x) {
            optimize1(x, M, I1, S1, I2, S2, l.iI, l.iS, logit.psi1.c[t, ], 
                      logit.psi2.c[t, ], alpha.c[t, ], equal.variance)
        }
                
        ## Make sure parameters fall within optimization bounds
        s1.u.old <- max(s.lower, s1.u.old)
        s2.u.old <- max(s.lower, s2.u.old)
        s.u.old <- max(s.lower, s.u.old)
                
        s1.c.old <- max(s.lower, s1.c.old)
        s2.c.old <- max(s.lower, s2.c.old)
        s.c.old <- max(s.lower, s.c.old)
                
        if (equal.variance == FALSE) {
            res1.u <- bobyqa(c(s1.u.old, s2.u.old, s.u.old, mu.u.old, delta.u.old), 
                             MLE1.u, lower = MLE1.lower.u, upper = MLE1.upper.u)
            res1.c <- bobyqa(c(s1.c.old, s2.c.old, s.c.old, mu.c.old, 0), MLE1.c, 
                             lower = MLE1.lower.c, upper = MLE1.upper.c)
                  
            s1.u[t] <- res1.u$par[1]
            s2.u[t] <- res1.u$par[2]
            s.u[t] <- res1.u$par[3]
            mu.u[t] <- res1.u$par[4]
            delta.u[t] <- res1.u$par[5]
            
            s1.c[t] <- res1.c$par[1]
            s2.c[t] <- res1.c$par[2]
            s.c[t] <- res1.c$par[3]
            mu.c[t] <- res1.c$par[4]
            delta.c[t] <- res1.c$par[5]
            
        } else {
            res1.u <- bobyqa(c(s1.u.old, s.u.old, mu.u.old, delta.u.old), MLE1.u, 
                             lower = MLE1.lower.u[c(1, 3, 4, 5)],
                             upper = MLE1.upper.u[c(1, 3, 4, 5)])
            res1.c <- bobyqa(c(s1.c.old, s.c.old, mu.c.old, 0), MLE1.c,
                             lower = MLE1.lower.c[c(1, 3, 4, 5)],
                             upper = MLE1.upper.c[c(1, 3, 4, 5)])
            
            s1.u[t] <- res1.u$par[1]
            s2.u[t] <- res1.u$par[1]
            s.u[t] <- res1.u$par[2]
            mu.u[t] <- res1.u$par[3]
            delta.u[t] <- res1.u$par[4]
            
            s1.c[t] <- res1.c$par[1]
            s2.c[t] <- res1.c$par[1]
            s.c[t] <- res1.c$par[2]
            mu.c[t] <- res1.c$par[3]
            delta.c[t] <- res1.c$par[4]
        }
                
                
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~#~ Step 3: Estimate logit psi's and alpha ~#~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        
        ## x[1] = logit_psi1 x[2] = logit_psi2 x[3] = alpha
        
        ## This is the third step in the optimization process where we estimate the MLEs
        ## of logit(psi1), logit(psi2), and alpha, based on the MLEs of delta, mu, sigma,
        ## sigma1, sigma2 computed in the previous stage.  See optimize2 for more details
        
        ## .u's always referred to unconstrained MLE .c's always referred to constrained
        ## MLE
                
        for (k in seq_len(M)) {
            MLE2.u <- function(x) {
                optimize2(x, k, I1, S1, I2, S2, l.iI, l.iS, delta.u[t], mu.u[t], 
                              s1.u[t], s2.u[t], s.u[t])
            }
                  
            MLE2.c <- function(x) {
                optimize2(x, k, I1, S1, I2, S2, l.iI, l.iS, delta.c[t], mu.c[t], 
                          s1.c[t], s2.c[t], s.c[t])
            }
                  
            res2.u <- optim(c(logit.psi1.u.old[k], logit.psi2.u.old[k], alpha.u.old[k]), 
                            MLE2.u)
            res2.c <- optim(c(logit.psi1.c.old[k], logit.psi2.c.old[k], alpha.c.old[k]), 
                            MLE2.c)
                  
            logit.psi1.u[t + 1, k] <- res2.u$par[1]
            logit.psi2.u[t + 1, k] <- res2.u$par[2]
            alpha.u[t + 1, k] <- res2.u$par[3]
            
            logit.psi1.c[t + 1, k] <- res2.c$par[1]
            logit.psi2.c[t + 1, k] <- res2.c$par[2]
            alpha.c[t + 1, k] <- res2.c$par[3]
            
        }
                
                
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~~~~~#~ Step 4: Evaluate log likelihood ~#~~~~~~~~~~~~~~~ ##
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        
        ll.new.u <- loglikelihood(M, I1, S1, I2, S2, l.iI, l.iS,
                                  logit.psi1.u[t, ], logit.psi2.u[t, ], alpha.u[t, ],
                                  s1.u[t], s2.u[t], s.u[t], mu.u[t], delta.u[t])
        ll.new.c <- loglikelihood(M, I1, S1, I2, S2, l.iI, l.iS,
                                  logit.psi1.c[t, ], logit.psi2.c[t, ], alpha.c[t, ],
                                  s1.c[t], s2.c[t], s.c[t], mu.c[t], delta.c[t])
                
        if ((abs(ll.new.u - ll.old.u) < tol) &
            (abs(ll.new.c - ll.old.c) < tol)) {
            break
        }
                
        ll.old.u <- ll.new.u
        ll.old.c <- ll.new.c
                
        logit.psi1.u.old <- logit.psi1.u[t + 1, ]
        logit.psi2.u.old <- logit.psi2.u[t + 1, ]
        alpha.u.old <- alpha.u[t + 1, ]
        
        logit.psi1.c.old <- logit.psi1.c[t + 1, ]
        logit.psi2.c.old <- logit.psi2.c[t + 1, ]
        alpha.c.old <- alpha.c[t + 1, ]
                
        s1.u.old <- s1.u[t]
        s2.u.old <- s2.u[t]
        s.u.old <- s.u[t]
        mu.u.old <- mu.u[t]
        delta.u.old <- delta.u[t]
                
        s1.c.old <- s1.c[t]
        s2.c.old <- s2.c[t]
        s.c.old <- s.c[t]
        mu.c.old <- mu.c[t]
        delta.c.old <- delta.c[t]
        
    }
            
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~#~ Step 5: Save all values ~#~~~~~~~~~~~~~~~ ##
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    
    ## Sometimes test.statistics can be negative due to numerical reasons.  In these
    ## cases, set them to a very small number, so they correspond to p-value = 1.
    test.stat <- -2 * (ll.new.c - ll.new.u)
    test.stat[test.stat < 0] <- 10^(-6)
    
    output <- list(test.stat, delta.u[t], mu.u[t], mu.c[t],
                   s1.u[t], s2.u[t], s.u[t], s1.c[t], s2.c[t], s.c[t],
                   sigmoid(logit.psi1.u[t + 1, ]), 
                   sigmoid(logit.psi2.u[t + 1, ]),
                   alpha.u[t + 1, ],
                   sigmoid(logit.psi1.c[t + 1, ]),
                   sigmoid(logit.psi2.c[t + 1, ]),
                   alpha.c[t + 1, ], M, t)
    ##exonList[iExon], as.vector(I1 - pseudocount), as.vector(S1 - pseudocount), 
    ##as.vector(I2 - pseudocount), as.vector(S2 - pseudocount))
    
    output
}

