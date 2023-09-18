#define ARMA_NO_DEBUG

#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec tmpU(double tmpV1,double tmpV2,double tmpV3,double tmpV4,double tmpV5,double tmpV6) {
    
    const int npow = 6;
    arma::vec tmpRes(npow);
    tmpRes.fill(0);
    
    tmpRes(0) = tmpV1;
    tmpRes(1) = tmpV1 * tmpV1 - tmpV2;
    tmpRes(2) = tmpV1 * tmpV1 * tmpV1 - 3 * tmpV2 * tmpV1 + 2 * tmpV3;
    tmpRes(3) = pow(tmpV1,4) - 6 * tmpV2 * pow(tmpV1,2) + 3 * pow(tmpV2,2) + 8 * tmpV3 * tmpV1 - 6 * tmpV4;
    tmpRes(4) = pow(tmpV1,5) - 10 * tmpV2 * pow(tmpV1,3) + 15 * pow(tmpV2,2) * tmpV1 + 20 * tmpV3 * pow(tmpV1,2) - 20 * tmpV3 * tmpV2 - 30 * tmpV4 * tmpV1 + 24 * tmpV5;
    tmpRes(5) = pow(tmpV1,6) - 15 * tmpV2 * pow(tmpV1,4) + 40 * tmpV3 * pow(tmpV1,3) + 45 * pow(tmpV1,2) * pow(tmpV2,2) -90 * pow(tmpV1,2) * tmpV4 - 120 * tmpV1 * tmpV2 * tmpV3 + 144 * tmpV1 * tmpV5 - 15 * pow(tmpV2,3) + 90 * tmpV2 * tmpV4 + 40 * pow(tmpV3,2) - 120 * tmpV6;
    
    return tmpRes;
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec computeU(arma::mat designX) {
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    const int q = p;// (p-1);
    const int npow = 6;
    arma::mat crossX(n,q);
    crossX.fill(0);
    
    arma::mat tmpRes(npow,q);
    tmpRes.fill(0);
    
    arma::mat tmpSum(p,1);
    tmpSum.fill(0);
    
    //caculate sample covariance based on X variables
    arma::mat sam_cov(p,p);
    sam_cov.fill(0);
    sam_cov = arma::cov(designX);
    
    arma::mat tmpUinf0(p,1);

    for(int i = 0; i < p; i++) {
                arma::mat crossXtmp = designX.col(i);
                double tmpV1 = accu(crossXtmp);
                double tmpV2 = accu(pow(crossXtmp,2));
                double tmpV3 = accu(pow(crossXtmp,3));
                double tmpV4 = accu(pow(crossXtmp,4));
                double tmpV5 = accu(pow(crossXtmp,5));
                double tmpV6 = accu(pow(crossXtmp,6));
                
                tmpRes.col(i) = tmpU(tmpV1,tmpV2,tmpV3,tmpV4,tmpV5,tmpV6);
                tmpSum.row(i) = pow(tmpV1,2)/n;

                //Calculate the statistics for SPU(inf)
                tmpUinf0.row(i) = tmpSum.row(i)/sam_cov(i,i);
    }
    
    arma::vec SPU(npow + 1);
    SPU.fill(0);
    
    arma::vec tmpSPU(npow);
    tmpSPU.fill(0);
    tmpSPU = arma::sum(tmpRes,1);
    
    for (int i = 0; i < npow; i++) {
        SPU(i) = tmpSPU(i);
    }
    
    // Standardize
    SPU(0) = SPU(0) / n;
    SPU(1) = SPU(1) / (n * (n-1));
    SPU(2) = SPU(2) / (n * (n-1) * (n-2) );
    SPU(3) = SPU(3) / (n * (n-1) * (n-2) * (n-3));
    SPU(4) = SPU(4) / (n * (n-1) * (n-2) * (n-3) * (n-4) );
    SPU(5) = SPU(5) / (n * (n-1) * (n-2) * (n-3) * (n-4) * (n-5));
    SPU(6) = tmpUinf0.max();
    
    return SPU;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat pvalSPU0(arma::mat XUs, arma::mat yresid, arma::mat powV, int nperm) {
    const int n = XUs.n_rows;
    const int p = XUs.n_cols;
    
    const int npow = powV.n_rows;
    // containers
    arma::mat T0s1(npow, nperm);
    T0s1.fill(0);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat r10 = shuffle(yresid);
        arma::mat U01(n,p);
        U01.fill(0);
        for (int k=0; k<p; k++){
            U01.col(k) =XUs.col(k) % r10; //column-wise product
        }
        
        T0s1.col(i)=computeU(U01);
    }
    //Rcpp::List res;
    //res["T0"] =T0s1;
    
    return(T0s1);
}

