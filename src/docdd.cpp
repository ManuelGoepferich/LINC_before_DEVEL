#include <Rcpp.h>
using namespace Rcpp;

// Compute the dice distance


// [[Rcpp::export]]

NumericMatrix docdd(NumericMatrix inmatrix,
                  NumericMatrix outmatrix, double th) {
  
  int coln = inmatrix.ncol();
  int rown = inmatrix.nrow();
  NumericVector x = inmatrix.nrow();
  NumericVector y = inmatrix.nrow();
  NumericVector fillshared = inmatrix.nrow();
  NumericVector fillsymetr = inmatrix.nrow();
  
  for(int nn = 0; nn < coln; ++nn) {
    for(int mm = 0; mm < rown; ++mm) {
      x[mm] = inmatrix(mm,nn);
    }
    for(int n = 0; n < coln ; ++n) {
      for(int m = 0; m < rown; ++m) {
        y[m] = inmatrix(m,n);
        fillshared[m] = fillsymetr[m] = 0;  
        
        if((x[m] > th) && (y[m] > th)){
          fillshared[m] = 1;
        } else if ((x[m] > th) ^ (y[m] > th)){
          fillsymetr[m] = 1;
        }
      }
      double shared = sum(fillshared);
      double symetr = sum(fillsymetr);
      double cddist = (symetr)/(2*shared + symetr);
      
      outmatrix(nn, n) = cddist;
    }
  }
  return outmatrix;
}