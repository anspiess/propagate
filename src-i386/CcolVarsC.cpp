#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP C_colVarsC(SEXP mat) {
  NumericMatrix x(mat);
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(ncol);
  
  /* first pass => mean */ 
    for (int i = 0; i < ncol; i++) {
      double n = 0;
      double sum1 = 0, sum2 = 0;
      double mean1;    
      
      for (int j = 0; j < nrow; j++) {
        double val = x(j, i);
        if (!ISNA(val)) {
          n++;
          sum1 += val;
        }
      }
      mean1 = sum1/n;      
      
      /* second pass => variance */
      for (int k = 0; k < nrow; k++) {
        double val = x(k, i);
        if (!ISNA(val)) sum2 = sum2 + (val - mean1) * (val - mean1);  
      }
      
      double var = sum2/(n - 1);        
      out[i] = var;
    }
  
  return (out);
}