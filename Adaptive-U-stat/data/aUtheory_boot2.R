apval_aU_boot <- function(Y, X,cov = NULL, pow = c(1:6,Inf), resample = "boot", model= "gaussian",n.perm = 5000){
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    if(is.null(cov)) {
        r <- Y - mean(Y)
        U <- t(X) %*% r
        U <- U / n
        X.new = sweep(X,MARGIN=1,r,`*`)
    } else {
        regression.data <- as.data.frame(cbind(trait = Y, cov))
        if(is.null(colnames(cov))) {
            colnames(regression.data) <- c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(regression.data) <- c("trait", colnames(cov))
        }
        
        fit <- glm(trait~.,family = model, data = regression.data)
        pis <- fitted.values(fit)
        rownames(regression.data) %in% names(pis)
        r <- Y - pis
        U <- t(X) %*% r
        U <- U / n
        X.new <- sweep(X,MARGIN=1,r,`*`)

    }
    
    len = length(pow)
    
    ## calculate the expecatation, variance, and covariance of L(gamma)
    parametric.boot <- matrix(NA, n.perm, len)
    if(resample == "boot") {
        for (b in 1:n.perm){
            set.seed(b)
            if (is.null(cov)) {
                Y0 <- sample(Y, length(Y))
                ##  Null score vector:
                r0 <- Y - mean(Y)
                X.new.tmp <- sweep(X,MARGIN=1,r0,`*`)
                U.tmp <- t(X) %*% r0
                U.tmp <- U.tmp / n
                parametric.boot[b,]<- SPU_stat2(X.new.tmp,U.tmp)
            } else {
                
                if( model == "gaussian" ) {
                    Y0 <- pis + sample(fit$residuals, n, replace = F )
                    tdat0 <- data.frame(trait=Y0, cov)
                    fit0 <-  glm(trait ~., data = tdat0)
                    yfits0 <- fitted.values(fit0)
                    r0 <- Y0 - yfits0
                    X.new.tmp <- sweep(X,MARGIN=1,r0,`*`)
                    U.tmp <- t(X) %*% r0
                    U.tmp <- U.tmp / n
                    parametric.boot[b,]<- SPU_stat2(X.new.tmp,U.tmp)
                } else {
                    ## with nuisance parameters:
                    Y0 <- rbinom(n,1, prob = pis)
                    # for(i in 1:n) Y0[i] <- sample(c(1,0), 1, prob=c(pis[i], 1-pis[i]) )
                    tdat0<-data.frame(trait=Y0, cov)
                    fit0<-glm(trait~., family=model, data=tdat0)
                    yfits0<-fitted.values(fit0)
                    r0 <- Y0 - yfits0
                    X.new.tmp <- sweep(X,MARGIN=1,r0,`*`)
                    U.tmp <- t(X) %*% r0
                    U.tmp <- U.tmp / n
                    parametric.boot[b,]<- t(computeU(X.new.tmp)) #SPU_stat2(X.new.tmp,U.tmp)
                }
            }
        }
    } else {
        for (b in 1:n.perm) {
            r0 <- sample(r, length(r))
            X.new.tmp <- sweep(X,MARGIN=1,r0,`*`)
            U.tmp <- t(X) %*% r0
            U.tmp <- U.tmp / n
            parametric.boot[b,]<- t(computeU(X.new.tmp)) #SPU_stat2(X.new.tmp,U.tmp)
        }
    }
    
    ##observed statistics
    T = t(computeU(X.new.tmp))   #SPU_stat2(X.new,U)
    
    L.e <- colMeans(parametric.boot)
    V <- diag(var(parametric.boot))

    L.e = L.e[1:6]
    V = V[1:6]
    T = T[1:6]
    T.final = T / sqrt(V)
    pval =  2 * (1 - pnorm(abs(T.final)))
    
    sam.cov = cov(X.new)
    diag.sam.cov <- diag(sam.cov)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
    L <- n * max(U^2/diag.sam.cov)
    stan.L.inf <- L- (2*log(p) - log(log(p)))
    pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(pi))
    
    pval = c(pval,pval.inf)

    pval.min = min(pval)
    aSPU = 1 - (1 - pval.min)^len
    
    pval.f = pchisq(  sum(- 2 * log(pval)) ,df = 2*7, lower.tail = FALSE)
    
    pvalue = c(pval, aSPU, pval.f )
    names(pvalue) = c("SPU1","SPU2","SPU3","SPU4","SPU5","SPU6","SPUInf","aUmin","aUf")
    
    Ts = c(T[1:6], L)
    
    pPerm0 = rep(NA,length(pow))
    T0s = parametric.boot
    tmp.f = rep(0,n.perm, 1)
    for ( j in 1:length(pow))
    {
        pPerm0[j] = round( sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm, digits = 8)
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm)
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
        tmp.f = tmp.f - 2*log(P0s)
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    
    T.Fisher = sum(-2*log(pPerm0+1/(n.perm))) #fisher statistic from permutation
    PaU_F <- (sum(T.Fisher <=  tmp.f) +1)/(n.perm+1)
    
    Ts <- c(Ts, min(pPerm0))
    
    pvs <- c(pPerm0, Paspu, PaU_F)
    names(Ts) <- c(paste("SPU", pow, sep=""),"aUmin")
    names(pvs) = c("SPU1","SPU2","SPU3","SPU4","SPU5","SPU6","SPUInf","aUmin","aUf")
    
    p.second= c(pval[1:6], pPerm0[7]) #asymptoc normal+permutation
    p.second.min= 1 - (1 -min(p.second) )^len 
    p.second.fisher= pchisq(  sum(- 2 * log(p.second)) ,df = 2*7, lower.tail = FALSE)
    
    out = list(Ts = Ts, pvalue1 = pvalue, pvalue3 = pvs, pvalues_second=c(p.second.min, p.second.fisher) )
    return(out)
}