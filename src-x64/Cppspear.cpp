#include <Rcpp.h>
using namespace Rcpp;

// Spearman complete
// [[Rcpp::export]]

NumericMatrix Cppspear(NumericMatrix xinmatrix,
                                    NumericMatrix yinmatrix) {
  
  int xnrow = xinmatrix.nrow();
  int ynrow = yinmatrix.nrow();
  int len   = xinmatrix.ncol();            
  NumericMatrix outmatrix(xnrow, ynrow);
  
  for(int k = 0; k  < xnrow ; ++k) {
    for(int j = 0; j < ynrow; ++j) {
      
      // to access the rows
      NumericVector xrow = xinmatrix.ncol();
      NumericVector yrow = yinmatrix.ncol();
      for(int h = 0; h < len; ++h) {
        xrow[h] = xinmatrix(k,h);
        yrow[h] = yinmatrix(j,h);
      }
      
      // remove NA values
      LogicalVector x_na = !is_na(xrow);
      LogicalVector y_na = !is_na(yrow);
      for(int l = 0; l < len; ++l){
        if(y_na[l] && x_na[l]){
          x_na[l] = TRUE;
        } else {
          x_na[l] = FALSE;
        }
      }
      
      // ranks for x value vector
      NumericVector x_clean = xrow[x_na];
      NumericVector x_sorted = clone(x_clean);
      NumericVector x_ranks = clone(x_clean);
      std::sort(x_sorted.begin(), x_sorted.end());
      int len_na   = x_clean.size();
      
      for(int i = 0; i < len_na; ++i){
        double e_sum = 0; double e_count = 0;
        for(int ii = 0; ii < len_na; ++ii) {
          if(x_clean[i] == x_sorted[ii]){
            e_sum = e_sum + (ii + 1);
            e_count = e_count + 1;
          }
        }
        x_ranks[i] = (e_sum/e_count);
      }   
      
      // ranks for y value vector
      NumericVector y_clean = yrow[x_na];
      NumericVector y_sorted = clone(y_clean);
      NumericVector y_ranks = clone(y_clean);
      std::sort(y_sorted.begin(), y_sorted.end());
      
      for(int i = 0; i < len_na; ++i){
        double e_sum = 0; double e_count = 0;
        for(int ii = 0; ii < len_na; ++ii) {
          if(y_clean[i] == y_sorted[ii]){
            e_sum = e_sum + (ii + 1);
            e_count = e_count + 1;
          }
        }
        y_ranks[i] = (e_sum/e_count);
      }   
      
      // diff rank
      double difxy = 0;
      for(int d = 0; d < len_na; ++d) {
        difxy += pow((x_ranks[d] - y_ranks[d]),2);
      }
      // final spearman computation
      double qspear;double nfraction;
      nfraction = len_na*(pow(len_na, 2)-1);
      outmatrix(k,j) = 1 - ((6*(difxy))/nfraction);
      
    }
  }
  return outmatrix;
}
  