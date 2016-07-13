#include <Rcpp.h>
using namespace Rcpp;

// OFFICIAL esd inner function
// do the t.test for esd

// [[Rcpp::export]]
NumericVector doesd(NumericVector query,
                                 NumericVector wquery, int rmax,
                                 NumericVector querypos,
                                 NumericVector rmaxvec,
                                 NumericVector mdiff,
                                 double alpha) {
  int n = query.size(); int nn = wquery.size();
  std::vector<double> rvalues (rmax,0);
  std::vector<double> pvalues (rmax,0);
  NumericVector wdiff = mdiff;
  
  // loop over i max outliers
  for(int i = 0; i < rmax ; ++i) {
    double meanq = mean(query);
    double sdq = sd(query);  
    
    for(int el = 0; el  < n ; ++el) {
      mdiff[el] = fabs((query[el] - meanq));
    }
    for(int ell = 0; ell  < nn ; ++ell) {
      wdiff[ell] = fabs((wquery[ell] - meanq));
      if(wdiff[ell] == max(mdiff)){
        double wpos = ell; querypos[i] = (wpos + 1);
      }
    }
    double rsingle = max(mdiff)/sdq;
    rvalues[i] = rsingle;
    int maxpos = which_max(mdiff);
    query.erase(maxpos); mdiff.erase(maxpos);
    
    // lambda
    int dof = (n - (i + 1) - 1); 
    double arange = (1 - (alpha/(2*(n - (i + 1) + 1))));
    n = n - 1;
    NumericVector quant = NumericVector::create(0, arange);
    NumericVector lambda = qt(quant, dof);
    pvalues[i] = lambda[1];
  }
  
  for(int rr = 0; rr < rmax; ++rr) {
    if(rvalues[rr] > pvalues[rr] ){
      rmaxvec[rr] = (rr + 1);
    }
  }
  // return candidate position  double maxnofr = max(rmaxvec);
  querypos.push_back(max(rmaxvec));
  return querypos;
}
