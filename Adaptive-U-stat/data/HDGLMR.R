HDGLMR <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), n.perm=1000){
    
    model = match.arg(model)
    #pow=c(2:8, Inf)
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r<-Y-mean(Y)
        U<-as.vector(t(X) %*% r)
    } else {
        tdat1 <- data.frame(trait=Y, cov)
        if(is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        fit1 <- glm(trait~.,family=model,data=tdat1)
        pis <- fitted.values(fit1)
        #r<- (Y - pis)/(pis * (1-pis))
        r<- Y - pis
        U<-t(X) %*% r
    }
    
    ##observed statistics
    M= 0
    for (i in 1:n) {
        M = M + r[i]^2 * t(as.matrix(X[i,])) %*% as.matrix(X[i,])
    }
    Ts =sum(U^2) - M
    
    TsC = HDGLMC(t(X),as.matrix(r),n.perm)
    pPerm0=  (sum(TsC$T0s1 >= Ts[1,1])+1) / (n.perm + 1)
    
    res = list(Ts = Ts, pvs = pPerm0)
    res

}


Goeman_test <- function(Y, X, cov = NULL, model=c("gaussian","binomial"), n.perm=1000){
    
    model = match.arg(model)
    #pow=c(2:8, Inf)
    n <- length(Y)
    if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
    k <- ncol(X)
    
    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r<-Y-mean(Y)
    } else {
        tdat1 <- data.frame(trait=Y, cov)
        if(is.null(colnames(cov))) {
            colnames(tdat1) = c("trait", paste("cov",1:dim(cov)[2],sep=""))
        } else {
            colnames(tdat1) = c("trait", colnames(cov))
        }
        fit1 <- glm(trait~.,family=model,data=tdat1)
        pis <- fitted.values(fit1)
        #r<- (Y - pis)/(pis * (1-pis))
        r<- Y - pis
    }
    
    ##observed statistics
    Ts = (t(r) %*% X %*% t(X) %*% r) /(t(r) %*% diag(X) %*% diag(X) %*% r)
    
    TsC = GeomanC(X,as.matrix(r),as.matrix(diag(X)),n.perm)
    pPerm0=  (sum(TsC$T0s1 >= Ts[1,1])+1) / (n.perm + 1)
    
    res = list(Ts = Ts, pvs = pPerm0)
    res
    
}
