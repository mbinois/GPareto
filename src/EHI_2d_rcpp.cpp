#include <Rcpp.h>
using namespace Rcpp;

// // [[Rcpp::export]]
// Hypervolume from : 
//  (c) Michael Emmerich and Andre Deutz, LIACS, Leiden University, 2010
//  emmerich@liacs.nl, deutz@liacs.nl
//
// title Hypervolume 2d
// param P Two lines matrix of the corresponding points
// param x Reference point
// [[Rcpp::export]]

double hvolume2d_Rcpp(NumericMatrix S, double x1, double x2) {
  double h = 0;
  
  // version initiale
  for(int i = 0; i < S.nrow(); i++){
      if(i==0){
        h += (x1 - S(i,0)) * (x2 - S(i,1));
      }
      else{
        h += (x1 - S(i,0)) *  (S(i-1,1) - S(i,1));
      }
    }

// int k = S.nrow();
//  for(int i = 0; i < S.nrow(); i++){
//      if(i == 0){
//        //std::cout << "i = " << i << " " << S(k - i - 1,0) << " " << S(k - i - 1,1) << std::endl;
//        h += (x1 - S(k - i - 1,0)) * (x2 - S(k - i - 1,1));
//      }
//      else{
//        //std::cout << "i = " << i << " " << S(k - i - 1,0) << " " << S(k - i,1) << " " <<  S(k - i - 1,1) << std::endl;
//        h += (x1 - S(k - i - 1,0)) *  (S(k - i,1) - S(k - i - 1,1));
//      }
//  }
    
  return(h);
}

// [[Rcpp::export]]
double exipsi_Rcpp(double a, double b, double m, double s){
  NumericVector tmp(1, (b-m)/s);
  NumericVector res(1);

  res = s*dnorm(tmp) + (a-m)*pnorm(tmp);
  
  return(res(0));
}


// [[Rcpp::export]]
NumericVector EHI_2d_wrap_Rcpp(NumericMatrix P,NumericVector r, NumericMatrix mu, NumericMatrix s){
  
 // determine all lower left corner cell coordinates
  
  
  int n = mu.nrow();
  int k = P.nrow();
  NumericVector resu(n);
  for (int candidat = 0; candidat < n; candidat++){  
  
  NumericVector c2(P(_,1));
  std::sort(c2.begin(), c2.end());
  NumericVector c1(P(_,0));
  
  NumericMatrix ccc(k+1,k+1);
  double fMax2 = 0;
  double fMax1 = 0;
  double cL1 = 0;
  double cL2 = 0;
  double cU1 = 0;
  double cU2 = 0;
  double sPlus = 0;
  double Psi1 = 0;
  double Psi2 = 0;
  double GaussCDF1 = 0;
  double GaussCDF2 = 0;
  
  for(int i = 0; i<=k; i++){
    for(int j = 0; j <= (k-i); j++){
      
      // c1(i), c2(j) are now the cell coordinates according Fig. 2   
      // For coordinate j determine hight fMax2 
      if (j==0) {
        fMax2=r(1);
      } else{
        fMax2 = c2(k-j);
      }
      
      // For coordinate i determine the width of the staircase fMax1 
      if (i==0) {
        fMax1 = r(0);
      } else {
        fMax1 = c1(k-i);
      }
    
      // get cell coordinates
      if (j==0) {
        cL1 = R_NegInf;
      } else{
        cL1 = c1(j-1);
      }
      if (i==0) {
        cL2= R_NegInf;
      } else {
        cL2 = c2(i-1);
      }
      if (j==k) {
        cU1 = r(0);
      } else{ 
        cU1 = c1(j);
      }
      if (i==k) {
        cU2 = r(1);
      } else { 
        cU2 = c2(i);
      }
      
      // SM = points that are dominated or equal to upper cell bound
      
      // find the size of SM
      int tmp = 0;
      for(int m = 0; m < k; m++){
        if(cU1 <= P(m,0) && cU2 <= P(m,1))
          tmp++;
      }
      
      if(tmp == 0){
        sPlus = 0;
      }
      else{
        NumericMatrix SM(tmp,2);
      
        int ind1 = 0;
        for (int m = 0; m < k; m++){
          if (cU1 <= P(m,0) && cU2 <= P(m,1)){
            SM(ind1, 0) = P(m,0);
            SM(ind1, 1) = P(m,1);
            ind1++;
          }
        }
        
        sPlus = hvolume2d_Rcpp(SM, fMax1,fMax2);
      }
      
      // Marginal integration over the length of a cell
      NumericVector candidateMu(mu(candidat,_));
      NumericVector candidateSigma(s(candidat,_));
      
      Psi1 = exipsi_Rcpp(fMax1,cU1,candidateMu(0),candidateSigma(0)) - exipsi_Rcpp(fMax1,cL1,candidateMu(0),candidateSigma(0));
      
      // Marginal integration over the height of a cell
      Psi2 = exipsi_Rcpp(fMax2,cU2,candidateMu[1],candidateSigma[1]) - exipsi_Rcpp(fMax2,cL2,candidateMu[1],candidateSigma[1]);
      
      // Cumulative Gaussian over length for correction constant
      NumericVector tmp1(1, (cU1-candidateMu(0))/candidateSigma(0));
      NumericVector tmp2(1, (cL1-candidateMu(0))/candidateSigma(0));
      
      tmp1 = pnorm(tmp1);
      tmp2 = pnorm(tmp2);
      
      GaussCDF1 = tmp1(0) - tmp2(0);
      
      // Cumulative Gaussian over length for correction constant
      tmp1 = (cU2-candidateMu[1])/candidateSigma[1];
      tmp1 = pnorm(tmp1);
      
      tmp2 = (cL2-candidateMu[1])/candidateSigma[1];
      tmp2 = pnorm(tmp2);
      
      GaussCDF2 = tmp1(0) - tmp2(0);
      
      // ExI Contribution from the actual cell
      ccc(i,j) = Psi1*Psi2-sPlus*GaussCDF1*GaussCDF2;
      
    }
  }
  
  double res = 0;
  for(int i = 0; i < k+1; i++){
    for(int j = 0; j < k+1;j++ ){
      res += std::max(ccc(i,j), 0.0);
    }
  }
  resu[candidat]= res;
}
  
  return(resu);
}

