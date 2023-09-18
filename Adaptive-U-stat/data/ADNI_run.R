rm(list =ls())
setwd("/home/panwei/wuxx0845/UTest/RealData/Pathway/ADNI1")
library(data.table)
#setwd("/Users/chong/Dropbox/aSPUtheory/RealData")

library(aSPU)
library(mvnfast)
library(truncnorm)
library(aSPU)
library(mvtnorm)
source("aSPUtheory.R")
library(HDGLM)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("HDGLMC.cpp")
source("HDGLMR.R")
source("aSPUtheory_boot.R")
source("aSPUtheory.R")
source("support.R")

args=(commandArgs(TRUE))
job = as.numeric(gsub("\\job=", "", args))

job = 15

pheno = read.csv(file="keyADNItables.csv", header=TRUE, sep=",")
pheno = pheno[,c("RID","PTID","ORIGPROT","DX.bl","AGE","PTGENDER","PTEDUCAT","ICV.bl","PTETHCAT","PTRACCAT","PTMARRY")]


pheno = pheno[!duplicated(pheno),]
pheno2 = pheno[pheno[,"ORIGPROT"]=="ADNI1",]
#pheno2 = pheno2[pheno2[,"PTETHCAT"]=="Not Hisp/Latino",]


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

GeneList = glist[glist$gene %in%pathway,]


GeneList[,1] = as.numeric(levels(GeneList[,1])[GeneList[,1]])
GeneList = GeneList[order(GeneList[,1]),]
GeneList[is.na(GeneList[,1]),1] = 23
index.set = unique(GeneList[,1])


chr = index.set[1]
raw = fread(paste("/home/panwei/wuxx0845/Imputed_ADNI1_ADNI2/ADNI1_processed/PathwayChr",chr,".raw",sep=""), header=T)
raw = as.data.frame(raw)

raw.id = raw[,2]


tmp = strsplit(raw.id,"_")
tmp = unlist(tmp)
tmp.len = length(tmp)
raw.id2 = paste(tmp[1:tmp.len %%4==2],tmp[1:tmp.len %%4==3],tmp[1:tmp.len %%4==0],sep="_")
rownames(raw) = raw.id2

raw[,2] = raw.id2
common.id = intersect(pheno2.id,raw.id2)


raw = raw[raw[,2] %in%common.id,]

raw = raw[,7:dim(raw)[2]]
gene = GeneList[GeneList[,1]==chr,]

snp.id = colnames(raw)
snp.id = gsub(".*:","",snp.id)
snp.id = gsub("\\_.*","",snp.id)

snp.id.num = as.numeric(snp.id)

snp.ind = matrix(NA,length(snp.id),dim(gene)[1])
for(index in 1:dim(gene)[1]) {
    snp.ind[,index] = (snp.id.num > (gene[index,2] - 20000) ) & (snp.id.num < (gene[index,3] + 20000))
}

snp.ind2 = rowSums(snp.ind,na.rm = T)
snp.ind2 = snp.ind2 > 0


final.data = raw[,snp.ind2]

index.set = index.set[-1]
index.set = index.set[index.set!=23] #We don't have X chromosome data

for(chr in index.set ) {
    
   raw = fread(paste("/home/panwei/wuxx0845/Imputed_ADNI1_ADNI2/ADNI1_processed/PathwayChr",chr,".raw",sep=""), header=T)
   raw = as.data.frame(raw)
   
   raw.id = raw[,2]
   
   
   tmp = strsplit(raw.id,"_")
   tmp = unlist(tmp)
   tmp.len = length(tmp)
   raw.id2 = paste(tmp[1:tmp.len %%4==2],tmp[1:tmp.len %%4==3],tmp[1:tmp.len %%4==0],sep="_")
   rownames(raw) = raw.id2
   
   raw[,2] = raw.id2
   common.id = intersect(pheno2.id,raw.id2)
   
   
   raw = raw[raw[,2] %in%common.id,]
   
   raw = raw[,7:dim(raw)[2]]
   gene = GeneList[GeneList[,1]==chr,]
   
   snp.id = colnames(raw)
   snp.id = gsub(".*:","",snp.id)
   snp.id = gsub("\\_.*","",snp.id)
   
   snp.id.num = as.numeric(snp.id)
   
   snp.ind = matrix(NA,length(snp.id),dim(gene)[1])
   for(index in 1:dim(gene)[1]) {
       snp.ind[,index] = (snp.id.num > (gene[index,2] - 20000) ) & (snp.id.num < (gene[index,3] + 20000))
   }
   
   snp.ind2 = rowSums(snp.ind,na.rm = T)
   snp.ind2 = snp.ind2 > 0
   
   
    geno = raw[,snp.ind2]

   
    final.data = merge(final.data, geno,by="row.names",all.x=TRUE)
    rownames(final.data) = final.data[,1]
    final.data= final.data[,-1]
}

final.data = final.data[,!duplicated(t(final.data))]

final.data = final.data[,colSums(is.na(final.data)) < 1]
save.name = paste("pathway",job,".rds",sep="")
saveRDS(final.data,save.name)

#Get covariates
pheno2$Y = as.numeric(pheno2[,"DX.bl"]=="CN")
pheno2$ICV.bl2 = pheno2$ICV.bl/10^6
pheno2$PTHAND = as.numeric(pheno2$PTHAND==1)
pheno2$PTGENDER.y = as.numeric(pheno2$PTGENDER.y==1)
pheno.used = pheno2[,c("Y","AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
rownames(pheno.used) = rownames(pheno2)

final.data1 = na.omit(final.data)

rm(final.data)

final.data1 = as.matrix(final.data1)

cov = pheno.used[,c("Y","AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
cov = cov[rownames(cov) %in% rownames(final.data1),]
cov = cov[rownames(final.data1),]
Y = cov[,"Y"]
cov = cov[,c("AGE","PTEDUCAT","PTHAND","PTGENDER.y","ICV.bl2")]
cov = as.matrix(cov)
X = final.data1
rm(final.data1)

X <- apply(X, 2, function(w) w-mean(w))



final.res = matrix(NA,9,7)

n.perm = 1000
pow = c(1:6,Inf)

# Store the number of genes and the number of SNPs
final.res[1,3] = dim(X)[1]
final.res[1,4] = dim(X)[2]
final.res[1,5] = dim(GeneList)[1]

###
pow = c(1:6,Inf)

final.res[2,] = aSPU(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1000)$pvs
if (min(final.res[2,]) < 5e-3){final.res[2,] = aSPU(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e4)$pvs}
if (min(final.res[2,]) < 5e-4){final.res[2,] = aSPU(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e5)$pvs}
if (min(final.res[2,]) < 5e-5){final.res[2,] = aSPU(Y, X, cov = cov, resample = "perm",model = "binomial", pow = pow, n.perm = 1e6)$pvs}



tmp = apval_aSPU_boot(Y,X,cov, pow = c(1:5,Inf), resample = "boot", model= "binomial",n.perm = 10000)
final.res[3,] = tmp$pvalue1
final.res[4,] = tmp$pvalue3

if (min(final.res[4,]) < 5e-4)
{
    tmp = apval_aSPU_boot(Y,X,cov, pow = c(1:5,Inf), resample = "boot", model= "binomial",n.perm = 100000)
    final.res[3,] = tmp$pvalue1
    final.res[4,] = tmp$pvalue3
}

save.name = paste("Res",job,".rds",sep="")
saveRDS(final.res,save.name)


tmp = apval_aSPU_boot(Y,X,cov, pow = c(1:5,Inf), resample = "boot", model= "binomial",n.perm = 10000)
final.res[8,] = tmp$pvalue1
final.res[9,] = tmp$pvalue3

if (min(final.res[9,]) < 5e-4)
{
    tmp = apval_aSPU_boot(Y,X,cov, pow = c(1:5,Inf), resample = "boot", model= "binomial",n.perm = 100000)
    final.res[8,] = tmp$pvalue1
    final.res[9,] = tmp$pvalue3
}

save.name = paste("Res",job,".rds",sep="")
saveRDS(final.res,save.name)
