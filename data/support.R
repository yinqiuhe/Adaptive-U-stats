
perm <- function(n,k){choose(n,k) * factorial(k)}

#calculate the U test statistc
SPU_stat <- function(U,S,cov.est) {
    n = dim(S)[1]
    p = dim(S)[2]
    len = 6
    out = matrix(NA,1,len)
    V.1 = colSums(S)
    V.2 = colSums(S^2)
    V.3 = colSums(S^3)
    V.4 = colSums(S^4)
    V.5 = colSums(S^5)
    
    U1 = sum(V.1) / perm(n,1)
    U2 = sum(V.1 * V.1 - V.2) / perm(n,2)
    U3 = sum(V.1 * V.1 * V.1 - 3 * V.2 * V.1 + 2 * V.3) / perm(n,3)
    U4 = sum(V.1 * V.1 * V.1 * V.1 - 6 * V.2 * V.1 * V.1 + 3 * V.2 * V.2 + 8 * V.3 * V.1 - 6 * V.4) / perm(n,4)
    U5 = sum(V.1 * V.1 * V.1 * V.1 * V.1 - 10 * V.2 * V.1 * V.1 * V.1 + 15 * V.2 * V.2 * V.1 + 20 * V.3 * V.1 * V.1 - 20 * V.3 * V.2 - 30 * V.4 * V.1 + 24 * V.5) / perm(n,5)
    
    var.1 = sum(cov.est) / choose(n,1)
    var.2 = sum(cov.est^2) / choose(n,2)
    var.3 = sum(cov.est^3) / choose(n,3)
    var.4 = sum(cov.est^4) / choose(n,4)
    var.5 = sum(cov.est^5) / choose(n,5)

    T = c(U1,U2,U3,U4,U5)
    V = c(var.1,var.2,var.3,var.4,var.5)
    T.final = T / sqrt(V)
    pval =  2 * (1 - pnorm(abs(T.final)))
    
    diag.sam.cov <- diag(cov.est)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
    L <- n * max(U^2/diag.sam.cov)
    stan.L.inf <- L- (2*log(p) - log(log(p)))
    pval.inf <- 1 - exp(-exp(-stan.L.inf/2)/sqrt(pi))
    
    pval = c(pval,pval.inf)
    pval.min = min(c(pval))
    aSPU = 1 - (1 - pval.min)^len
    
    pvalue = c(pval, aSPU)
    names(pvalue) = c("SPU1","SPU2","SPU3","SPU4","SPU5","SPUInf","aSPU")
    names(T.final) = c("SPU1","SPU2","SPU3","SPU4","SPU5")
    
    out = list(pval = pvalue, T = T.final)
    return(out)
}


SPU_stat2 <- function(S,U) {
    n = dim(S)[1]
    len = 5
    out = matrix(NA,1,len)
    V.1 = colSums(S)
    V.2 = colSums(S^2)
    V.3 = colSums(S^3)
    V.4 = colSums(S^4)
    V.5 = colSums(S^5)
    
    U1 = sum(V.1) / perm(n,1)
    U2 = sum(V.1 * V.1 - V.2) / perm(n,2)
    U3 = sum(V.1 * V.1 * V.1 - 3 * V.2 * V.1 + 2 * V.3) / perm(n,3)
    U4 = sum(V.1 * V.1 * V.1 * V.1 - 6 * V.2 * V.1 * V.1 + 3 * V.2 * V.2 + 8 * V.3 * V.1 - 6 * V.4) / perm(n,4)
    U5 = sum(V.1 * V.1 * V.1 * V.1 * V.1 - 10 * V.2 * V.1 * V.1 * V.1 + 15 * V.2 * V.2 * V.1 + 20 * V.3 * V.1 * V.1 - 20 * V.3 * V.2 - 30 * V.4 * V.1 + 24 * V.5) / perm(n,5)
    
    cov.est <- cov(S)
    diag.sam.cov <- diag(cov.est)
    diag.sam.cov[diag.sam.cov <= 10^(-10)] <- 10^(-10)
    Uinf <- n * max(U^2/diag.sam.cov)
    T = c(U1,U2,U3,U4,U5, Uinf)
    
    names(T) = c("SPU1","SPU2","SPU3","SPU4","SPU5","SPUInf")
    
    return(T)
}


best.band <- function(sam, bandwidth, cv.fold= 5, norm= "F"){
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        sam.train.cov <- cov(sam.train)
        sam.test.cov <- cov(sam.test)
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


best.band2 <- function(sam, bandwidth, cv.fold= 5, norm= "F"){
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        sam.train.cov <- t(sam.train) %*% sam.train / (dim(sam.train)[1]-1)
        
        sam.test.cov <- t(sam.test) %*% sam.test / (dim(sam.test)[1]-1)
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


best.band3 <- function(sam, cov.new, bandwidth, cv.fold= 5, norm= "F"){
    p <- dim(sam)[2]
    n <- dim(sam)[1]
    fold.size <- round(n/cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    for(i in 1:cv.fold){
        if(i == cv.fold){
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):n]
        }else{
            temp.idx <- sam.idx[((i - 1)*fold.size + 1):(i*fold.size)]
        }
        sam.train <- sam[-temp.idx,]
        sam.test <- sam[temp.idx,]
        cov.train <- cov.new[-temp.idx,]
        cov.test <- cov.new[temp.idx,]
        
        V11 <- t(cov.train) %*% cov.train / (dim(cov.train)[1]-1)
        V12 <- t(cov.train) %*% sam.train / (dim(cov.train)[1]-1)
        V22 <- t(sam.train)%*% sam.train / (dim(sam.train)[1]-1)
        sam.train.cov <- V22 - t(V12) %*% solve(V11) %*% V12
        
        
        V11 <- t(cov.test) %*% cov.test / (dim(cov.test)[1]-1)
        V12 <- t(cov.test) %*% sam.test / (dim(cov.test)[1]-1)
        V22 <- t(sam.test)%*% sam.test / (dim(sam.test)[1]-1)
        sam.test.cov <- V22 - t(V12) %*% solve(V11) %*% V12
        
        for(j in 1:n.bandwidth){
            sam.train.cov.band <- sam.train.cov
            sam.train.cov.band[abs(row(sam.train.cov.band) - col(sam.train.cov.band)) > bandwidth[j]] <- 0
            diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
        }
    }
    diff.norm <- colMeans(diff.norm)
    best <- which(diff.norm == min(diff.norm))
    best <- best[1]
    return(bandwidth[best])
}


SPU_E <- function(gamma, n, cov){
    if(gamma %% 2 == 1){
        e <- 0
    }
    if(gamma %% 2 == 0){
        ga.half <- gamma/2
        e <- (factorial(gamma)*sum(diag(cov)^ga.half))/(2^ga.half * factorial(ga.half) * n^{ga.half})
    }
    return(e)
}


SPU_E2 <- function(gamma, n, cov){
    out = 0
    if(gamma %% 2 == 0) {
        ga.half <- gamma/2
        e <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half * factorial(ga.half) * n^{ga.half})
        out = e%*%e
    }
    return(out)
}

SPU_E3 <- function(s,t, n, cov){
    out = 0
    p <- dim(cov)[1]
    gamma = s
    ga.half <- gamma/2
    e1 <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half*factorial(ga.half) * n^{ga.half})
    gamma = t
    ga.half <- gamma/2
    e2 <- (factorial(gamma)*diag(cov)^ga.half)/(2^ga.half*factorial(ga.half) * n^{ga.half})
    out = e1 %*% e2
    
    return(out)
}


c_space <- function(s, t){
    c1.vec <- c2.vec <- c3.vec <- NULL
    for(c3 in 0:min(c(t, s))){
        for(c1 in 0:floor((t - c3)/2)) {
            for(c2 in 0:floor((s - c3)/2)) {
                if( t - 2*c1 - c3 == 0 & s - 2*c2 - c3  == 0){
                    c1.vec <- c(c1.vec, c1)
                    c2.vec <- c(c2.vec, c2)
                    c3.vec <- c(c3.vec, c3)
                }
            }
        }
    }
    return(data.frame(c1 = c1.vec, c2 = c2.vec, c3 = c3.vec))
}



SPU_Var <- function(gamma,n, cov){
    p <- dim(cov)[1]
    if(gamma == 1){
        var <- (1/n)*(t(rep(1, p)) %*% cov %*% rep(1, p))
        var <- as.numeric(var)
    } else if(gamma %%2 ==0) {
        P1 <- SPU_E(2*gamma, n, cov)
        P2 <- -SPU_E2(gamma, n, cov)
        C <- c_space(gamma, gamma)
        C = C[C[,3] >0,]
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- (factorial(gamma))^2*sum(mat)
            D <- (n^gamma * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P2 + P3
    } else {
        P1 = SPU_E(2* gamma, n, cov)
        C <- c_space(gamma, gamma)
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- (factorial(gamma))^2*sum(mat)
            D <- (n^gamma * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P3
    }
    return(var)
}



SPU_Cov <- function(s,t,n, cov){
    p <- dim(cov)[1]
    var = 0
    if(s %% 2 ==0) {
        P1 <- SPU_E(s+t, n, cov)
        P2 <- -SPU_E3(s,t, n, cov)
        C <- c_space(s, t)
        C = C[C[,3] >0,]
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- factorial(s) * factorial(t) * sum(mat)
            D <- (n^((s+t)/2) * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P2 + P3
    } else {
        P1 <- SPU_E(s+t, n, cov)
        C <- c_space(s, t)
        n.case <- dim(C)[1]
        P3 <- 0
        
        diags <- diag(cov)
        p <- length(diags)
        mat1 <- matrix(rep(diags, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diags, p), p, p, byrow = TRUE)
        for(i in 1:n.case){
            c1 <- C$c1[i]
            c2 <- C$c2[i]
            c3 <- C$c3[i]
            
            mat <- mat1^(c1)*mat2^(c2)*cov^(c3)
            diag(mat) <- 0
            N <- factorial(s) * factorial(t) * sum(mat)
            D <- (n^((s+t)/2) * factorial(c1)*factorial(c2)*factorial(c3)*2^(c1 + c2))
            P3 <- P3 + N/D
        }
        var <- P1 + P3
    }
    return(var)
}





