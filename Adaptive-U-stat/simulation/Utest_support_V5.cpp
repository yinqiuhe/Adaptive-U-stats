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
double testU2(arma::mat designX, int i, int j){
//function that computes the full U2 statistic

    const int n = designX.n_rows;
    //const int p = designX.n_cols;
    
    //arma::vec singleU2(1);
    double singleU2;
    
    //compute the single U(2) given i,j
    double smj1 = accu(designX.col(i));
    double smj2 = accu(designX.col(j));
    double smj1j1 = accu(pow(designX.col(i),2));
    double smj2j2 = accu(pow(designX.col(j),2));
    double smj1j2 = accu(designX.col(i) % designX.col(j));
    double smj1j1j2 = accu( pow(designX.col(i),2) %  designX.col(j));
    double smj1j2j2 = accu( designX.col(i) %  pow(designX.col(j),2) );
    double smj1j1j2j2 =  accu( pow(designX.col(i),2) % pow( designX.col(j) ,2) );
    
    double s3term = smj1j2*(smj1*smj2-smj1j2)-( smj1j1j2*smj2 - smj1j1j2j2 ) - ( smj1j2j2*smj1 - smj1j1j2j2 );
    double s4term = (smj1*smj1 - smj1j1 )*( smj2*smj2 - smj2j2 ) -2*( smj1j2*smj1j2 - smj1j1j2j2 ) - 4*s3term;

    //singleU2(0) =  (smj1j2*smj1j2 -smj1j1j2j2)/ (n * (n-1))-2 * s3term/(n*(n-1)*(n-2)) + s4term/(n*(n-1)*(n-2)*(n-3));
    singleU2 =  (smj1j2*smj1j2 -smj1j1j2j2)/ (n * (n-1))-2 * s3term/(n*(n-1)*(n-2)) + s4term/(n*(n-1)*(n-2)*(n-3));
    
    return singleU2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat centerdata(arma::mat designX){
    const int n = designX.n_rows;
    arma::mat colmean = mean(designX , 0);
    arma::vec allonev(n);
    allonev.ones();
    arma::mat centerm = kron(colmean,allonev);
    arma::mat centeredx = designX -centerm;
    return centeredx;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec computeU(arma::mat designX) {
//compute the U-statistics with U2 using full computation and the others with centered data
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    arma::mat designXnon = designX; //non-centered X
    designX = centerdata(designX); //center X data
    
    const int q = p * (p-1)/2;// (p-1);
    const int npow = 6;
    arma::mat crossX(n,q);
    crossX.fill(0);
    
    //arma::vec U2(1);
    double U2;
    
    arma::mat tmpRes(npow,q);
    tmpRes.fill(0);
    
    int tmpl = 0;
    for(int i = 0; i < p; i++) {
        for (int j = 0; j < i; j++) {
            //if (j != i) {
            //compute the general U statistics
                arma::mat crossXtmp = designX.col(i) % designX.col(j);
                double tmpV1 = accu(crossXtmp);
                double tmpV2 = accu(pow(crossXtmp,2));
                double tmpV3 = accu(pow(crossXtmp,3));
                double tmpV4 = accu(pow(crossXtmp,4));
                double tmpV5 = accu(pow(crossXtmp,5));
                double tmpV6 = accu(pow(crossXtmp,6));
            
                tmpRes.col(tmpl) = tmpU(tmpV1,tmpV2,tmpV3,tmpV4,tmpV5,tmpV6);

                tmpl++;
           
            //compute SPU(2) separately
                U2 = U2+ testU2(designXnon, i ,j );
            
           // }
        }
    }
    
    arma::vec SPU(npow + 1);
    SPU.fill(0);
    
    arma::vec tmpSPU(npow);
    tmpSPU.fill(0);
    tmpSPU = arma::sum(tmpRes,1);
    
    for (int i = 0; i < npow; i++) {
        SPU(i) = 2 * tmpSPU(i);
    }
    //compute SPU(2) separately
    SPU(1) = 2* U2;
    
    // Standardize
    SPU(0) = SPU(0) / n;
    //SPU(1) = SPU(1) / (n * (n-1));
    SPU(2) = SPU(2) / (n * (n-1) * (n-2) );
    SPU(3) = SPU(3) / (n * (n-1) * (n-2) * (n-3));
    SPU(4) = SPU(4) / (n * (n-1) * (n-2) * (n-3) * (n-4) );
    SPU(5) = SPU(5) / (n * (n-1) * (n-2) * (n-3) * (n-4) * (n-5));
    
    //Calculate the statistics for SPU(inf)
    arma::mat correlation(p,p);
    correlation.fill(0);
    
    correlation = arma::cor(designX,designX);
    correlation = arma::abs(correlation);
    
    arma::mat corDiag(p,p);
    corDiag.eye();
    
    correlation = correlation - corDiag;
    SPU(6) = arma::max(arma::max(correlation));

    return SPU;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec computeUold(arma::mat designX) {
//the old U-statistics computation without any centering of data
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    const int q = p * (p-1)/2;// (p-1);
    const int npow = 6;
    arma::mat crossX(n,q);
    crossX.fill(0);
    
    arma::mat tmpRes(npow,q);
    tmpRes.fill(0);
    
    int tmpl = 0;
    for(int i = 0; i < p; i++) {
        for (int j = 0; j < i; j++) {
            //if (j != i) {
            arma::mat crossXtmp = designX.col(i) % designX.col(j);
            double tmpV1 = accu(crossXtmp);
            double tmpV2 = accu(pow(crossXtmp,2));
            double tmpV3 = accu(pow(crossXtmp,3));
            double tmpV4 = accu(pow(crossXtmp,4));
            double tmpV5 = accu(pow(crossXtmp,5));
            double tmpV6 = accu(pow(crossXtmp,6));
            
            tmpRes.col(tmpl) = tmpU(tmpV1,tmpV2,tmpV3,tmpV4,tmpV5,tmpV6);
            
            tmpl++;
            
            // }
        }
    }
    
    arma::vec SPU(npow + 1);
    SPU.fill(0);
    
    arma::vec tmpSPU(npow);
    tmpSPU.fill(0);
    tmpSPU = arma::sum(tmpRes,1);
    
    for (int i = 0; i < npow; i++) {
        SPU(i) = 2 * tmpSPU(i);
    }
    
    // Standardize
    SPU(0) = SPU(0) / n;
    SPU(1) = SPU(1) / (n * (n-1));
    SPU(2) = SPU(2) / (n * (n-1) * (n-2) );
    SPU(3) = SPU(3) / (n * (n-1) * (n-2) * (n-3));
    SPU(4) = SPU(4) / (n * (n-1) * (n-2) * (n-3) * (n-4) );
    SPU(5) = SPU(5) / (n * (n-1) * (n-2) * (n-3) * (n-4) * (n-5));
    
    //Calculate the statistics for SPU(inf)
    arma::mat correlation(p,p);
    correlation.fill(0);
    
    correlation = arma::cor(designX,designX);
    correlation = arma::abs(correlation);
    
    arma::mat corDiag(p,p);
    corDiag.eye();
    
    correlation = correlation - corDiag;
    SPU(6) = arma::max(arma::max(correlation));
    
    return SPU;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List pvalSPU(arma::mat designX, int nperm) {
// compute p-values of all U-statistics using permutation of data
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    arma::vec SPU = computeU(designX);
    const int npow = 7;
    arma::mat T0s1(npow,nperm);
    T0s1.fill(0);
    for (int k = 0; k< nperm; k++) {
        arma::mat shuffledX(n,p);
        shuffledX.fill(0);
        
        for (int i = 0; i< p; i++) {
            shuffledX.col(i) = shuffle(designX.col(i));
            //  T0s1.col(i) =shuffledX.col(0);
        }
        
        T0s1.col(k) =computeU(shuffledX);
    }
    
    Rcpp::List res;
    res["T0"] =T0s1;
    res["Test"] = SPU;
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
double testInf(arma::mat designX) {
// compute U(inf)
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    double tInf;
    
    //Calculate the statistics for SPU(inf)
    arma::mat correlation(p,p);
    correlation.fill(0);
    
    correlation = arma::cor(designX,designX);
    correlation = arma::abs(correlation);
    
    arma::mat corDiag(p,p);
    corDiag.eye();
    
    correlation = correlation - corDiag;
    tInf = arma::max(arma::max(correlation));
    
    tInf = n*pow(tInf,2);
    
    return tInf;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List pvalInf(arma::mat designX, int nperm) {
// compute p-value of U(inf) using permutation
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    double tInf0 = testInf(designX);

    arma::vec T0inf(nperm);
    T0inf.fill(0);
    for (int k = 0; k < nperm; k++) {
        arma::mat shuffledX(n,p);
        shuffledX.fill(0);
        
        for (int i = 0; i< p; i++) {
            shuffledX.col(i) = shuffle(designX.col(i));
            //  T0s1.col(i) =shuffledX.col(0);
        }

        T0inf(k) = testInf(shuffledX);
    }
    
    Rcpp::List res;
    res["T0"] =T0inf;
    res["Test"] = tInf0;
    return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
double testInf2(arma::mat designX){
// compute U(inf) with standardization on denominator
    
    //const int p = designX.n_cols;
    const int n = designX.n_rows;
    arma::mat cen_X = centerdata(designX);
    arma::mat sigma_hat = cen_X.t() * cen_X/n;
    arma::mat cent_X_sq = cen_X % cen_X;
    arma::mat shat = cent_X_sq.t() * cent_X_sq / n -sigma_hat % sigma_hat;
    arma::mat t_pert = abs( sigma_hat/ sqrt(shat/n) );
    t_pert.diag().zeros();
    
    double t_inf = arma::max(arma::max(t_pert));
    
    return t_inf;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List pvalInf2(arma::mat designX, int nperm) {
// compute p-value of U(inf) with standardization on denominator by permutation
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    
    double tInf0 = testInf2(designX);
    
    arma::vec T0inf(nperm);
    T0inf.fill(0);
    for (int k = 0; k < nperm; k++) {
        arma::mat shuffledX(n,p);
        shuffledX.fill(0);
        
        for (int i = 0; i< p; i++) {
            shuffledX.col(i) = shuffle(designX.col(i));
            //  T0s1.col(i) =shuffledX.col(0);
        }
        
        T0inf(k) = testInf2(shuffledX);
    }
    
    Rcpp::List res;
    res["T0"] =T0inf;
    res["Test"] = tInf0;
    return res;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec compute_sd(arma::mat designX){
// compute the SDs of U-statistics with p^4 order summation
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    int gamma_num = 6;
    
    arma::vec var_est(gamma_num);
    arma::vec tmpures(gamma_num);
    
    var_est.fill(0);
    tmpures.fill(0);
    
    
    for (int j1= 0; j1< (p-1); j1++) {
        for (int j2= j1 + 1 ; j2 < p; j2++) {
            for (int j3 = 0; j3<(p-1); j3++) {
                for (int j4 = j3 + 1; j4 <p; j4++) {
                    arma::mat crossXtmp = designX.col(j1) % designX.col(j2) % designX.col(j3) % designX.col(j4);
                    double tmpV1 = accu(crossXtmp);
                    double tmpV2 = accu(pow(crossXtmp,2));
                    double tmpV3 = accu(pow(crossXtmp,3));
                    double tmpV4 = accu(pow(crossXtmp,4));
                    double tmpV5 = accu(pow(crossXtmp,5));
                    double tmpV6 = accu(pow(crossXtmp,6));
                    
                    tmpures = tmpU(tmpV1,tmpV2,tmpV3,tmpV4,tmpV5,tmpV6);
                    
                    var_est(0) = var_est(0) + tmpures(0);
                    var_est(1) = var_est(1) + tmpures(1);
                    var_est(2) = var_est(2) + tmpures(2);
                    var_est(3) = var_est(3) + tmpures(3);
                    var_est(4) = var_est(4) + tmpures(4);
                    var_est(5) = var_est(5) + tmpures(5);
                }
            }
        }
    }

    var_est(0) = sqrt( 4 * var_est(0)/(n * n));
    var_est(1) = sqrt( 4 * 2 * var_est(1) / pow(n * (n-1),2) );
    var_est(2) = sqrt( 4 * 6 * var_est(2) / pow(n * (n-1) * (n-2),2) );
    var_est(3) = sqrt( 4 * 24 * var_est(3) / pow(n * (n-1) * (n-2) * (n-3),2) );
    var_est(4) = sqrt( 4 * 120 * var_est(4) / pow(n * (n-1) * (n-2) * (n - 3) * (n - 4),2) );
    var_est(5) = sqrt( 4 * 720 * var_est(5) / pow(n * (n-1) * (n-2) * (n - 3) * (n - 4) * (n - 5),2) );

    return(var_est);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec compute_sd_ind(arma::mat designX){
// compute the SDs of U-statistics with p^2 order summation
    const int n = designX.n_rows;
    const int p = designX.n_cols;
    int gamma_num = 6;
    
    arma::vec var_est(gamma_num);
    arma::vec tmpures(gamma_num);
    
    var_est.fill(0);
    tmpures.fill(0);
    
    
    for (int j1= 0; j1< (p-1); j1++) {
        for (int j2= j1 + 1 ; j2 < p; j2++) {
            arma::mat crossXtmp = pow( designX.col(j1) % designX.col(j2), 2 ) ;
            double tmpV1 = accu(crossXtmp);
            double tmpV2 = accu(pow(crossXtmp,2));
            double tmpV3 = accu(pow(crossXtmp,3));
            double tmpV4 = accu(pow(crossXtmp,4));
            double tmpV5 = accu(pow(crossXtmp,5));
            double tmpV6 = accu(pow(crossXtmp,6));
            
            tmpures = tmpU(tmpV1,tmpV2,tmpV3,tmpV4,tmpV5,tmpV6);
            
            var_est(0) = var_est(0) + tmpures(0);
            var_est(1) = var_est(1) + tmpures(1);
            var_est(2) = var_est(2) + tmpures(2);
            var_est(3) = var_est(3) + tmpures(3);
            var_est(4) = var_est(4) + tmpures(4);
            var_est(5) = var_est(5) + tmpures(5);
        }
    }
    
    var_est(0) = sqrt( 4 * var_est(0)/(n * n));
    var_est(1) = sqrt( 4 * 2 * var_est(1) / pow(n * (n-1),2) );
    var_est(2) = sqrt( 4 * 6 * var_est(2) / pow(n * (n-1) * (n-2),2) );
    var_est(3) = sqrt( 4 * 24 * var_est(3) / pow(n * (n-1) * (n-2) * (n-3),2) );
    var_est(4) = sqrt( 4 * 120 * var_est(4) / pow(n * (n-1) * (n-2) * (n - 3) * (n - 4),2) );
    var_est(5) = sqrt( 4 * 720 * var_est(5) / pow(n * (n-1) * (n-2) * (n - 3) * (n - 4) * (n - 5),2) );
    
    return(var_est);
}
