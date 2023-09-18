#use permutation to compute p-value of U(inf), just rho
inf.test<-function(design.X, nperm=1000) 
{
  res = pvalInf(design.X,nperm)
  
  t.inf =  res$T0
  t.inf = t(t.inf)
  
  t0.inf = res$Test
  
  pval.inf = sum( abs(t0.inf) <= abs(t.inf) ) / nperm
  
  return(pval.inf)
}

#use permutation to compute p-value of U(inf), standardized rho
inf.test2<-function(design.X, nperm) 
{
  res = pvalInf2(design.X,nperm)
  t.inf =  res$T0
  t.inf = t(t.inf)
  t0.inf = res$Test
  pval.inf = sum( abs(t0.inf) <= abs(t.inf) ) / nperm
  return(pval.inf)
}


#function on computing 6 p-values of finite Us by using the variance estimator method
pcalc.Utest<-function(design_X){
  Ts=computeU(design_X) #U2 computes the complete test statistic, other Us computed by centered data
  design_X =  design_X- matrix(rep(colMeans(design_X),each = n),n,p) #center the same data
  SDs=compute_sd_ind(design_X) 
  tstat = Ts[1:6]/SDs 
  Upvals=pnorm(abs(tstat),lower.tail = FALSE)*2  
  return(Upvals)
}

#function on different combinations of p-values
#return 7 p-values of single Ustatistics; 6 different combination of p-values (min and Fisher for two different max-p-value or do not use max-p-value)
p.val.comb<-function(Upvals,max.pval1,max.pval2){
  pval = c(Upvals,max.pval1)
  aSPU = 1- (1-min(pval))^(length(pval))
  T.aSPU2 = sum(- 2 * log(pval))
  aSPU2 = pchisq(T.aSPU2,df = 2*7, lower.tail = FALSE)
  
  pval.tmp = pval[1:6]
  aSPU0 = 1- (1-min(pval.tmp))^(length(pval.tmp))
  T.aSPUF.finite = sum(- 2 * log(pval.tmp))
  aSPUF.f = pchisq(T.aSPUF.finite,df = 2*6, lower.tail = FALSE)
  
  pvals = c(Upvals,max.pval2)
  aSPUs = 1- (1-min(pvals))^(length(pvals))
  T.aSPU2s = sum(- 2 * log(pvals))
  aSPU2s = pchisq(T.aSPU2s,df = 2*7, lower.tail = FALSE)
  
  res=c( Upvals,max.pval1,max.pval2, aSPU,aSPU2,  aSPUs,aSPU2s,  aSPU0,aSPUF.f)
  names(res)=c("vSPU(1)","vSPU(2)","vSPU(3)","vSPU(4)","vSPU(5)","vSPU(6)","vSPU(Inf)1","vSPU(Inf)2","v.aSPU","v.aSPU_Fisher", "v.aSPU2","v.aSPU_Fisher2","v.aSPU_finite","v.aSPU_Fisher_finite")
  return(res)  
}













#other tested methods or functions



#function on computing p-values by using the variance estimator method
#max.pval: is the permutation p-value of maxU, if NULL, we compute maxU p-value by asymptotic in this function
#return: 6 p-values of finite U statistics, p-value of infinity U, min aggregation of finite p-values, Fisher aggregation of finite p-values, min aggregation of all p-values, Fisher aggregation of all p-values, 
pcalc.Utest.old.2<-function(design_X, max.pval=NULL){
  Ts=computeU(design_X) #U2 computes the complete test statistic, other Us computed by centered data
  design_X =  design_X- matrix(rep(colMeans(design_X),each = n),n,p) #center the same data
  SDs=compute_sd_ind(design_X) 
  tstat = Ts[1:6]/SDs 
  Upvals=pnorm(abs(tstat),lower.tail = FALSE)*2  
  
  if (is.null(max.pval)) {
    maxT <- n * Ts[7]^2
    stan.L.inf <- maxT - (4*log(p) - log(log(p)))
    max.pval <- 1 - exp(-exp(-stan.L.inf/2) / sqrt(8 * pi))
  }
  
  pval = c(Upvals,max.pval)
  
  aSPU = 1- (1-min(pval))^(length(pval))
 
  T.aSPU2 = sum(- 2 * log(pval))
  aSPU2 = pchisq(T.aSPU2,df = 2*7, lower.tail = FALSE)
  
  pval.tmp = pval[1:6]
  aSPU0 = 1- (1-min(pval.tmp))^(length(pval.tmp))
  
  T.aSPUF.finite = sum(- 2 * log(pval.tmp))
  aSPUF.f = pchisq(T.aSPUF.finite,df = 2*6, lower.tail = FALSE)
  
  res=c(Upvals,max.pval,aSPU0,aSPUF.f,aSPU,aSPU2)
  names(res)=c("vSPU(1)","vSPU(2)","vSPU(3)","vSPU(4)","vSPU(5)","vSPU(6)","vSPU(Inf)", "v.aSPU_finite","v.aSPU_Fisher_finite","v.aSPU","v.aSPU_Fisher")
  return(res)
} 


pcalc.Utest.old<-function(design_X, max.pval=NULL){
  Ts=computeUold(design_X) #use the old method without any centering on data
  SDs=compute_sd_ind(design_X)  #do not center the same data to compute variance
  tstat = Ts[1:6]/SDs 
  Upvals=pnorm(abs(tstat),lower.tail = FALSE)*2  
  
  if (is.null(max.pval)) {
    maxT <- n * Ts[7]^2
    stan.L.inf <- maxT - (4*log(p) - log(log(p)))
    max.pval <- 1 - exp(-exp(-stan.L.inf/2) / sqrt(8 * pi))
  }
  
  pval = c(Upvals,max.pval)
  
  aSPU = 1- (1-min(pval))^(length(pval))
  
  T.aSPU2 = sum(- 2 * log(pval))
  aSPU2 = pchisq(T.aSPU2,df = 2*7, lower.tail = FALSE)
  
  pval.tmp = pval[1:6]
  aSPU0 = 1- (1-min(pval.tmp))^(length(pval.tmp))
  
  T.aSPUF.finite = sum(- 2 * log(pval.tmp))
  aSPUF.f = pchisq(T.aSPUF.finite,df = 2*6, lower.tail = FALSE)
  
  res=c(Upvals,max.pval,aSPU0,aSPUF.f,aSPU,aSPU2)
  names(res)=c("vSPU(1)","vSPU(2)","vSPU(3)","vSPU(4)","vSPU(5)","vSPU(6)","vSPU(Inf)", "v.aSPU_finite","v.aSPU_Fisher_finite","v.aSPU","v.aSPU_Fisher")
  return(res)
} 



#multiplier bootstrap adapted from two-sample setting
pertInf.boot<-function(design.X, nperm=1000) 
{
  n=dim(design.X)[1]
  p=dim(design.X)[2]
  
  center.X= design.X- matrix(rep(colMeans(design.X),each = n),n,p)
  sigma.hat=t(center.X)%*%center.X/n
  cen.X.sq=center.X^2
  shat =t(cen.X.sq)%*%cen.X.sq/n- sigma.hat*sigma.hat 
  t.pert =abs( sigma.hat/sqrt(shat/n) )
  diag(t.pert)=0
  t.pert=max(t.pert)
  
  pert.res = matrix(0,1,nperm)
  for (ind in 1:nperm)
  {
    pert.res[ind] = perInf.test(center.X,sigma.hat,shat,n)
  }
  pval.inf=sum(t.pert<=pert.res)/nperm
  return(pval.inf)
}
perInf.test<-function(center.X, sigma.hat, shat, n)
{
  gn = rnorm(n)
  sigma.pert=  t(center.X)%*%diag(gn)%*%center.X/n-mean(gn)*sigma.hat
  t.pert =abs( sigma.pert/sqrt(shat/n) )
  diag(t.pert)=0
  t.pert=max(t.pert)
  return(t.pert)
}
asym.test<-function(design_X) #asymptotic rsult in Jiang
{
	  maxT <- testInf(design_X)
    stan.L.inf <- maxT - (4*log(p) - log(log(p)))
    max.pval <- 1 - exp(-exp(-stan.L.inf/2) / sqrt(8 * pi))
    return(max.pval)
}
asym.test2<-function(design_X) #asymptotic result in Liu tilde
{
	  maxT <- testInf(design_X)
    stan.L.inf <- maxT - (4*log(p) - log(log(p)))
    pex <- pchisq( ((4*log(p) - log(log(p)))+stan.L.inf),  df=1,  lower.tail=FALSE)
    max.pval <- 1 - exp(-(p^2-p)*pex/2)
    return(max.pval)
}
asym.test3<-function(design_X) #asymptotic result in Liu W
{
	n=dim(design_X)[1]
  p=dim(design_X)[2]

	nhf = floor(n/2)
	dX1 = design_X[1: nhf,]
	dX2 = design_X[(nhf+1):n,]
	ctX1 = dX1 - matrix(rep(colMeans(design_X),each = nhf),nhf,p)
	ctX2 = dX2 - matrix(rep(colMeans(design_X),each = (n-nhf)), (n-nhf), p )
  cX = design_X - matrix(rep(colMeans(design_X),each = n),n,p)
	   
  A = t(ctX1)%*%ctX1
  B = t(ctX2)%*%ctX2

	D = t(t(colSums(cX^2)))%*%colSums(cX^2)

	Lsq = (2*A^2+2*B^2)/D
	diag(Lsq)=0
	
	maxT <-  n*max(max(  Lsq ))-4*log(p)
    
  max.pval <- 1 - exp(-exp(-maxT/2)/2)
  return(max.pval)
}


