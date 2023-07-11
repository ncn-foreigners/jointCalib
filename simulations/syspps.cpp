#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector syspps_cpp(NumericVector x, int n) {
  int N = x.size();
  IntegerVector U = sample(N, N) - 1; // Subtract 1 because C++ indexing starts from 0
  NumericVector xx(N);
  NumericVector z(N);
  double x_sum = sum(x);
  double cumsum = 0;
  
  for (int i = 0; i < N; i++) xx[i] = x[U[i]];
  for (int i = 0; i < N; i++) {
    cumsum += xx[i];
    z[i] = n * cumsum / x_sum;
  }
  
  double r = unif_rand();
  IntegerVector s(n);
  int j = 0;
  
  for (int i = 0; i < N; i++) {
    while (z[i] >= r && j < n) {
      s[j] = U[i] + 1; // Adding 1 because we need to return 1-based index
      r += 1;
      j += 1;
    }
  }
  
  return s.sort();
}


