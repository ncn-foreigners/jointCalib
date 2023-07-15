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

// [[Rcpp::export]]
NumericVector inclusionprobabilities_cpp(NumericVector a, double n) {
  int len = a.size();
  double a_sum = sum(a);
  NumericVector pik1 = n * a / a_sum;
  NumericVector pik(len);
  std::vector<int> list1_index;
  
  for (int i = 0; i < len; ++i) {
    if (pik1[i] > 0) {
      list1_index.push_back(i);
      pik[i] = pik1[i];
    }
  }
  
  std::vector<int> list_index;
  for (int i : list1_index) {
    if (pik[i] >= 1) {
      list_index.push_back(i);
    }
  }
  
  int l = list_index.size();
  int l1;
  
  if (l > 0) {
    l1 = 0;
    while (l != l1) {
      double x_sum = 0.0;
      for (int i = 0; i < len; ++i) {
        if (std::find(list_index.begin(), list_index.end(), i) == list_index.end()) {
          x_sum += pik[i];
        }
      }
      
      for (int i = 0; i < len; ++i) {
        if (std::find(list_index.begin(), list_index.end(), i) == list_index.end()) {
          pik[i] = (n - l) * pik[i] / x_sum;
        } else {
          pik[i] = 1;
        }
      }
      
      l1 = l;
      list_index.clear();
      for (int i : list1_index) {
        if (pik[i] >= 1) {
          list_index.push_back(i);
        }
      }
      l = list_index.size();
    }
    for (int i : list1_index) {
      pik1[i] = pik[i];
    }
  }
  return pik1;
}

// [[Rcpp::export]]
NumericVector UPtille_cpp(NumericVector pik, double eps = 1e-6) {
  if (any(is_na(pik))) stop("there are missing values in the pik vector");
  
  int n = sum(pik);
  LogicalVector list = (pik > eps) & (pik < 1 - eps);
  NumericVector pikb = pik[list];
  int N = pikb.size();
  
  if (N < 1) stop("the pik vector has all elements outside of the range [eps,1-eps]");
  
  NumericVector sb(N, 1.0);
  NumericVector b(N, 1.0);
  
  for (int i = 0; i < (N - n); ++i) {
    NumericVector a = inclusionprobabilities_cpp(pikb, N - i);
    NumericVector v = 1 - a/b;
    b = a;
    NumericVector p = v * sb;
    p = NumericVector(cumsum(p));
    double u = R::runif(0, 1);
    int j = 0;
    while (u >= p[j] && j < N) ++j;
    sb[j] = 0;
  }
  
  pik[list] = sb;
  return pik;
}



