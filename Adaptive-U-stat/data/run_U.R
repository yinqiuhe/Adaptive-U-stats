rm(list =ls())
library(data.table)
# library(aSPU)
library(mvnfast)
library(truncnorm)
library(mvtnorm)
source("aSPUtheory.R")
library(HDGLM)
library(globaltest)
library('plink2R')

library(Rcpp)
library(RcppArmadillo)
sourceCpp("HDGLMC.cpp")
source("HDGLMR.R")
source("aSPUtheory_boot.R")
source("aSPU_perm.R")
sourceCpp("Utest_GLM.cpp")


source("support.R") #Chong's version of calculating U-statistics 1-5 and Inf
# source("aUtheory_boot.R") #Chong's version of bootstraping and permutation 1-5 and Inf
source("aUtheory_boot2.R") #my codes using 1-6 and Inf


args=(commandArgs(TRUE))
job = as.numeric(gsub("\\job=", "", args))



#####################################################################################
#read data and extract pathway

pheno = read.csv(file="keyADNItables.csv", header=TRUE, sep=",")
pheno = pheno[,c("RID","PTID","ORIGPROT","DX.bl","AGE","PTGENDER","PTEDUCAT","ICV.bl","PTETHCAT","PTRACCAT","PTMARRY")]

pheno = pheno[!duplicated(pheno),]
pheno2 = pheno[pheno[,"ORIGPROT"]=="ADNI1",]


pheno.handed = read.csv(file="PTDEMOG.csv", header=TRUE, sep=",")
pheno.handed = pheno.handed[pheno.handed[,"Phase"]=="ADNI1",]
pheno.handed = pheno.handed[,c("RID","SITEID","PTHAND","PTGENDER")] #,
pheno.handed = pheno.handed[!duplicated(pheno.handed),]
pheno.handed = pheno.handed[pheno.handed[,1] %in% pheno2$RID,]
pheno.handed = pheno.handed[pheno.handed[,"PTHAND"] != -4,]
pheno2 = merge(pheno2,pheno.handed,by = "RID")

pheno2 = pheno2[rowSums(is.na(pheno2))==0,]

pheno2.id = as.character(pheno2[,"PTID"])
rownames(pheno2) =pheno2.id

load("KEGGpathway.RData")

tmp = as.data.frame(lapply(KEGG.pathway, length))
tmp = tmp[,tmp>=10 & tmp<=2000] # we only analyze the pathway with 10~500 genes.
pathway.name = colnames(tmp)

pathway = KEGG.pathway[pathway.name[job]]

pathway = pathway[[1]]

glist = read.table("glist-hg19.txt", header=F)
glist$gene = as.character(glist[,4])


glist = read.table("glist-hg19.txt", header=F)
glist[,4] = as.character(glist[,4])



GeneList = glist[glist[,4] %in%pathway,]

GeneList[,1] = as.numeric(levels(GeneList[,1])[GeneList[,1]])
GeneList = GeneList[order(GeneList[,1]),]

GeneList = GeneList[!is.na(GeneList[,1]),]

geneset =GeneList
geneset[,2] = geneset[,2] - 20 *1000
geneset[,3] = geneset[,3] + 20 *1000


genotype = NULL



allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}



for(i in 1:dim(geneset)[1]) {
    
    system(paste0("./plink --bfile ","/panfs/roc/groups/1/panwei/shared/11Imputed_ADNI1_ADNI2_1000G_MichiganServer_By_Zhiyuan/ADNI1/chr",geneset[i,1],
    " --chr ",geneset[i,1]," --from-kb ",geneset[i,2]/1e3," --to-kb ",
    geneset[i,3]/1e3," --maf 0.05 --hwe 0.05 --make-bed --out ", geneset[i,4]))
    
    genos = tryCatch(
    {read_plink(paste0(geneset[i,4]),impute="avg")}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    if(!is.null(genos)) {
        
        qc = allele.qc(  genos$bim[,5]  ,genos$bim[,6] , genos$bim[,5] , genos$bim[,6] )
        
        if ( sum(!qc$keep) > 0 ) {
            genos$bim = genos$bim[qc$keep,]
            genos$bed = genos$bed[,qc$keep]
        }
        
        tmp.geno = genos$bed #353
        
        genotype = cbind(genotype,tmp.geno)
        system(paste0("rm ",geneset[i,4],"*"))
    }
}
#### extract X0:the un-weighted SNPs and Y
### remove ambiguous SNPs

final.data = genotype
final.data = final.data[,!duplicated(t(final.data))]

final.data = final.data[,colSums(is.na(final.data))<50]

raw.id = rownames(final.data)


tmp = strsplit(raw.id,"_")
tmp = unlist(tmp)
tmp.len = length(tmp)
raw.id2 = paste(tmp[1:tmp.len %%4==2],tmp[1:tmp.len %%4==3],tmp[1:tmp.len %%4==0],sep="_")
rownames(final.data) = raw.id2

pheno3 = pheno2
#Get covariates
pheno2 = pheno2[pheno2[,"DX.bl"]!="AD", ] #=="CN"|pheno2[,"DX.bl"]=
pheno2$Y = as.numeric(pheno2[,"DX.bl"]=="CN")
pheno2$ICV.bl2 = pheno2$ICV.bl/10^6
pheno2$PTGENDER.y = as.numeric(pheno2$PTGENDER.x=="Male")
pheno.used = pheno2[,c("Y","AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
rownames(pheno.used) = rownames(pheno2)

final.data1 = na.omit(final.data)

final.data1 = as.matrix(final.data1)

cov = pheno.used[,c("Y","AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
cov = cov[rownames(cov) %in% rownames(final.data1),]

final.data1 = final.data1[rownames(final.data1) %in% rownames(cov),]

final.data1 = final.data1[rownames(cov),]
Y = cov[,"Y"]
cov = cov[,c("AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
colnames(cov) = c("AGE","PTEDUCAT","PTHAND","PTGENDER","ICV")

#E.gender = cov[,"PTGENDER.y"]
#E.gender[E.gender==0] = -1
#cov = cov[,c("AGE","PTEDUCAT","ICV.bl2")]


cov = as.matrix(cov)
#cov = cbind(cov,E.gender)
X = final.data1

#X <- apply(X, 2, function(w) w-mean(w))




#####################################################################################
#run the codes
final.res = matrix(NA,13,9)

n.perm = 1000
pow = c(1:6,Inf)

final.res[1,1:8] =  apval_aSPU_boot(Y, X, cov = cov, pow = pow,model = "binomial",  n.perm = 10000)$pvs

X.tmp = cbind(1,cov,X)
result=HDGLM_test(Y,X.tmp,nuisance=c(1:6),model="logistic")
final.res[2,1] = result$test_pvalue
hdglm.perm = HDGLMR(Y, X, cov = cov,model = "binomial", n.perm = 1000)$pvs
if(hdglm.perm<5e-3) {
    hdglm.perm = HDGLMR(Y, X, cov = cov,model = "binomial", n.perm = 1e4)$pvs
}

if(hdglm.perm<5e-4) {
    hdglm.perm = HDGLMR(Y, X, cov = cov,model = "binomial", n.perm = 1e5)$pvs
}
# 
# if(hdglm.perm<5e-5) {
#     hdglm.perm = HDGLMR(Y, X, cov = cov,model = "binomial", n.perm = 1e6)$pvs
# }
final.res[2,2] = hdglm.perm


# Store the number of genes and the number of SNPs
final.res[2,3] = dim(X)[1]
final.res[2,4] = dim(X)[2]
final.res[2,5] = dim(GeneList)[1]


# Goeman globaltest
colnames(X) = paste("X",1:dim(X)[2],sep = "")

dat.tmp = cbind(Y,cov,X)

alterative = as.formula(paste("Y ~ AGE + PTEDUCAT + PTHAND + PTGENDER + ICV+",paste(colnames(X),collapse ="+")))
res.global = gt(Y ~ AGE + PTEDUCAT + PTHAND + PTGENDER + ICV, alterative,data = dat.tmp, model = "logistic")
final.res[3,1] = p.value(res.global)


n.perm = 1000
pow = c(1:6,Inf)

final.res[4,1:8] = apval_aSPU(Y,X, cov = cov,pow = c(1:6, Inf), bandwidth = 10,model = "binomial",cov.est.type = "Method1")$pvalue
final.res[5,1:8] = apval_aSPU(Y,X, cov = cov,pow = c(1:6, Inf), bandwidth = 20,model = "binomial",cov.est.type = "Method1")$pvalue
final.res[7,1:8] = apval_aSPU(Y,X, cov = cov,pow = c(1:6, Inf), bandwidth = 30,model = "binomial",cov.est.type = "Method1")$pvalue


final.res[6,1:8] = aSPU_perm(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1000)$pvs
if (min(final.res[6,]) < 5e-3){final.res[6,] = aSPU_perm(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e4)$pvs}
if (min(final.res[6,]) < 5e-4){final.res[6,] = aSPU_perm(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e5)$pvs}
# if (min(final.res[6,]) < 5e-5){final.res[6,] = aSPU_perm(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e6)$pvs}



# tmp = apval_aU_boot(Y,X,cov, pow = c(1:6,Inf), resample = "boot", model= "binomial",n.perm = 1e3)
# final.res[8,1:9] = tmp$pvalue1 #pvalues from asymptotics 
# final.res[9,1:9] = tmp$pvalue3 #pvalues from parametric bootstrap
# final.res[10,1:2] = tmp$pvalues_second
# if (min(tmp$pvalue3) < 5e-3){
#     tmp = apval_aU_boot(Y,X,cov, pow = c(1:6,Inf), resample = "boot", model= "binomial",n.perm = 1e5)
#     final.res[9,1:9] = tmp$pvalue3 #pvalues from parametric bootstrap
#     final.res[10,3:4] = tmp$pvalues_second
# }

tmp = apval_aU_boot(Y,X,cov, pow = c(1:6,Inf), resample = "perm", model= "binomial",n.perm = 1e3)
final.res[11,1:9] = tmp$pvalue1 #pvalues from asymptotics 
final.res[12,1:9] = tmp$pvalue3 #pvalues from permutation
final.res[13,1:2] = tmp$pvalues_second
if (min(tmp$pvalue3) < 5e-3){
    tmp = apval_aU_boot(Y,X,cov, pow = c(1:6,Inf), resample = "perm", model= "binomial",n.perm = 1e5)
    final.res[12,1:9] = tmp$pvalue3 #pvalues from parametric bootstrap
    final.res[13 ,3:4] = tmp$pvalues_second
}

save.name = paste("rev_Res",job,".rds",sep="")
saveRDS(final.res,save.name)




