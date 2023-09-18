LWtest<-function(design_X)
{
  n<-dim(design_X)[1]
  p<-dim(design_X)[2]
  design_X  = design_X- matrix(rep(colMeans(design_X),each = n),n,p) #center the data
  S=t(design_X)%*%design_X/(n-1)
  #compute W
  trS=sum(diag(S))
  trSI= sum(sum( (S-diag(p))^2 )) #trace of matrix square is the summation of square of each element
  W=trSI/p-p*( trS/p)^2/(n-1) +p/(n-1)
  tstat=((n-1)*W-p-1)/2
  pval=pnorm(abs(tstat),lower.tail = FALSE)*2
  return(pval)
}


schott.method<-function(design_X)
{
  n<-dim(design_X)[1]
  p<-dim(design_X)[2]
  design_X  = design_X- matrix(rep(colMeans(design_X),each = n),n,p) #center the data
  samvar=matrix(0,p,1) #compute sample variances
  for(j in 1:p){ samvar[j]=sum(design_X[,j]^2)/(n-1) }
  design_X=design_X/matrix(rep(sqrt(samvar),each=n),n,p) #standardize data by sample SD
  
  tstat=0
  for(j1 in 2:p){
    for (j2 in 1:(j1-1)){
      tstat=tstat+(sum(design_X[,j1]*design_X[,j2])/(n-1))^2 #summation of squared correlation
  }}   
  vart=p*(p-1)*(n-1-1)/((n-1)^2*(n+2-1))
  tstat=(tstat-p*(p-1)/(2*(n-1)))/sqrt(vart)
  pval=pnorm(abs(tstat),lower.tail = FALSE)*2
  return(pval)
}