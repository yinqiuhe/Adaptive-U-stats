########################################################################################### 
## The function for testing \Sigma=I_p (Identity test)
## INPUT: 
##        X: is a n by p matrix, each row is a p-dim sample
##        alpha: significant level
## OUTPUT: 
##        NewStat: the test statistic value for testing identity from Chen, Zhang and Zhong (2010);
##        New: reject indicator for the test proposed by Chen, Zhang and Zhong (2010);  
########################################################################################### 

equality<-function(X,alpha)
{
n<-dim(X)[1]
p<-dim(X)[2]
## CZZ Test statistics
XXn<-X%*%t(X)
Y1<-sum(diag(t(X)%*%X))/n
Y3<-(sum(XXn)-Y1*n)/(n*(n-1))
Y2<-(sum(XXn^2)-sum(diag(XXn^2)))/(n*(n-1))
XXn2<-XXn%*%XXn
diagXX<-rep(diag(XXn),each=n-1)
offdiagXX<-as.logical(lower.tri(XXn)+upper.tri(XXn))
VecOffdiagXX<-matrix(XXn,n^2,1)[offdiagXX]
Y4<-(sum(XXn2)-sum(diag(XXn2))-2*sum(diagXX*VecOffdiagXX))/(n*(n-1)*(n-2))
Y5<-((n*(n-1)*Y3)^2-2*n*(n-1)*Y2-4*n*(n-1)*(n-2)*Y4)/(n*(n-1)*(n-2)*(n-3))
Tn1<-Y1-Y3
Tn2<-Y2-2*Y4+Y5
VCZ<-n*(Tn2/p-2*Tn1/p+1)/2
rej2<-0
if (VCZ>qnorm(1-alpha,0,1))
 {rej2<-1}
return(list(NewStat=VCZ,New=rej2,pval = pnorm(VCZ,lower.tail =FALSE) ))
}


###########################################################################################
## The function for testing \Sigma=s^2I_p (Sphericity test)
## INPUT: 
##        X: is a n by p matrix, each row is a p-dim sample
##        alpha: significant level
## OUTPUT: 
##        NewStat: the test statistic value for testing Sphericity from Chen, Zhang and Zhong (2010);
##        New: reject indicator for the Sphericity test proposed by Chen, Zhang and Zhong (2010);  
########################################################################################### 

sphericity<-function(X,alpha)
{
n<-dim(X)[1]
p<-dim(X)[2]
## New Test statistics
XXn<-X%*%t(X)
Y1<-sum(diag(t(X)%*%X))/n
Y3<-(sum(XXn)-Y1*n)/(n*(n-1))
Y2<-(sum(XXn^2)-sum(diag(XXn^2)))/(n*(n-1))
XXn2<-XXn%*%XXn
diagXX<-rep(diag(XXn),each=n-1)
offdiagXX<-as.logical(lower.tri(XXn)+upper.tri(XXn))
VecOffdiagXX<-matrix(XXn,n^2,1)[offdiagXX]
Y4<-(sum(XXn2)-sum(diag(XXn2))-2*sum(diagXX*VecOffdiagXX))/(n*(n-1)*(n-2))
Y5<-((n*(n-1)*Y3)^2-2*n*(n-1)*Y2-4*n*(n-1)*(n-2)*Y4)/(n*(n-1)*(n-2)*(n-3))
Tn1<-Y1-Y3
Tn2<-Y2-2*Y4+Y5
UCZ<-n*(p*Tn2/(Tn1^2)-1)/2
rej2<-0
if (UCZ>qnorm(1-alpha,0,1))
 {rej2<-1}
return(list(NewStat=UCZ,New=rej2,pval = pnorm(UCZ,lower.tail =FALSE)))
}




