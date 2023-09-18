source("CZZ_test.R")
source("Other_test.R")
source("var_est_V4.R")
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(BH)
sourceCpp("Utest_support_V5.cpp")

args=(commandArgs(TRUE))
job = as.numeric(gsub("\\job=", "", args))

outd=paste("/home//simulation/sim_p_1000_a2/sim1_p1000a2")
system(paste("mkdir -p ", outd, sep=""))

setwd("/home//simulation/sim_p_1000_a2/sim1_p1000a2")

#setting
#n=100, p=1000, |J_A|=10, rho=0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 (10 in total)

#########################################################################
#generate data functions, where JAs represents the size of alternative
library(MASS)
gen_indx<-function(p,JAs) #generated random posistions: first generate JAs positions from all p(p-1)/2 possibilities and then adjust the numbers
{
  allsamk=sort(sample(p*(p-1)/2,JAs))
  newsamall=allsamk
  newsamall
  for (k in 1:(p-1))
  {
    indx=which( (  (2*p-k)*(k-1)/2+1  )  <=  allsamk  & allsamk <= ( (2*p-k-1)*k/2 ) )
    newsamall[ indx ] =  newsamall[ indx ] + k*(k+1)/2
  }
  return(newsamall)
}
sim_data <- function(n,p,rho,JAs = 10/2) { #JAs is the number of nonzero entries in half matrix
  sigma = matrix(0,p,p)
  allinx=gen_indx(p,JAs)
  sigma[allinx]=rho
  sigma=sigma+t(sigma)
  diag(sigma) = 1
  mu = rep(0,p)
  X = mvrnorm(n,mu,sigma)
  
  return(X)
}
#########################################################################

set.seed(job)

n=100
p=1000
mu0=0

#final.res = matrix(NA,11,15)
final.res = matrix(NA,11,18)

######record computing time############
design_X = sim_data(n,p,0)

time = proc.time()
res.inf.1 = inf.test(design_X,nperm = 1000) #compute p-value of U(inf)
time2 = proc.time() - time
final.res[11,1] = time2[1]

time = proc.time()
res.inf.2 = inf.test2(design_X,nperm = 1000) #compute p-value of U(inf) with standardization denominator
time2 = proc.time() - time
final.res[11,2] = time2[1]

time = proc.time()
res.finite = pcalc.Utest(design_X) #compute p-value of finite Ustats
time2 = proc.time() - time
final.res[11,3] = time2[1]

time = proc.time()
chen.equal = equality(design_X,0.05)$pval
time2 = proc.time() - time
final.res[11,4] = time2[1]

time = proc.time()
chen.sphericity = sphericity(design_X,0.05)$pval
time2 = proc.time() - time
final.res[11,5] = time2[1]

time = proc.time()
lw.res = LWtest(design_X)
time2 = proc.time() - time
final.res[11,6] = time2[1]

time = proc.time()
schott.res = schott.method(design_X)
time2 = proc.time() - time
final.res[11,7] = time2[1]

res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) #return 8 p-values of single Ustatistics (6+1+1) and 6 different combinations
final.res[1,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

#######################################

design_X = sim_data(n,p,0.1)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[2,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.2)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[3,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.3)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[4,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.4)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[5,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.5)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[6,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.6)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[7,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.7)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[8,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.8)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[9,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.9)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[10,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

save.name = paste0("out_",job,".rds")
saveRDS(final.res,save.name)


#########################################################################
#second job

job = job+530
set.seed(job)


n=100
p=1000
mu0=0

final.res = matrix(NA,11,18)

######record computing time############
design_X = sim_data(n,p,0)

time = proc.time()
res.inf.1 = inf.test(design_X,nperm = 1000) #compute p-value of U(inf)
time2 = proc.time() - time
final.res[11,1] = time2[1]

time = proc.time()
res.inf.2 = inf.test2(design_X,nperm = 1000) #compute p-value of U(inf) with standardization denominator
time2 = proc.time() - time
final.res[11,2] = time2[1]

time = proc.time()
res.finite = pcalc.Utest(design_X) #compute p-value of finite Ustats
time2 = proc.time() - time
final.res[11,3] = time2[1]

time = proc.time()
chen.equal = equality(design_X,0.05)$pval
time2 = proc.time() - time
final.res[11,4] = time2[1]

time = proc.time()
chen.sphericity = sphericity(design_X,0.05)$pval
time2 = proc.time() - time
final.res[11,5] = time2[1]

time = proc.time()
lw.res = LWtest(design_X)
time2 = proc.time() - time
final.res[11,6] = time2[1]

time = proc.time()
schott.res = schott.method(design_X)
time2 = proc.time() - time
final.res[11,7] = time2[1]

res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) #return 8 p-values of single Ustatistics (6+1+1) and 6 different combinations
final.res[1,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

#######################################

design_X = sim_data(n,p,0.1)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[2,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.2)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[3,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.3)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[4,] =  c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.4)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[5,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.5)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[6,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)


design_X = sim_data(n,p,0.6)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[7,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.7)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[8,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.8)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[9,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

design_X = sim_data(n,p,0.9)
res.inf.1 = inf.test(design_X,nperm = 1000)
res.inf.2 = inf.test2(design_X,nperm = 1000)
res.finite = pcalc.Utest(design_X)
res1 = p.val.comb(res.finite,res.inf.1,res.inf.2) 
chen.equal = equality(design_X,0.05)$pval
chen.sphericity = sphericity(design_X,0.05)$pval
lw.res = LWtest(design_X)
schott.res = schott.method(design_X)
final.res[10,] = c(res1, chen.equal, chen.sphericity, lw.res, schott.res)

save.name = paste0("out_",job,".rds")
saveRDS(final.res,save.name)

