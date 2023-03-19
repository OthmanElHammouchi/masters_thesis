#ifndef HELPERS
#define HELPERS

#include <stdlib.h>
#include <Rcpp.h>

inline double** array_create(int n, int m) {
    double* data = (double*)malloc(n * m * sizeof(double));
    double** array = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) {
      array[i] = &(data[m * i]);
    }
    return(array);
}

inline double** array_from_rcpp(Rcpp::NumericMatrix array) {
  int n = array.rows();
  int m = array.cols();
  double** res = array_create(n, m);
  for (int i=0; i < n; i++) {
    for (int j=0; j < m; j++) {
      res[i][j] = array(i, j);
    }
  }
  return(res); 
}

inline double** array_transpose(int n, int m, double** array) {
  double** res = array_create(m, n);
  for (int i=0; i < n; i++) {
    for (int j=0; j < m; j++) {
      res[j][i] = array[i][j];
    }
  }
  return(res);
}

#endif