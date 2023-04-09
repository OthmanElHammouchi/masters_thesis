#ifndef HELPERS
#define HELPERS

#include <Rcpp.h>
#include <stdlib.h>

inline Rcpp::NumericMatrix na_to_zero(Rcpp::NumericMatrix triangle) {
  int n_dev = triangle.rows();
  for (int i = 0; i<n_dev; i++) {
    for (int j = 0; j<n_dev; j++) {
      if (triangle(i, j) == Rcpp::NA) triangle(i, j) = 0;
    }
  }
  return(triangle);
}

inline bool contains_str(Rcpp::CharacterVector vec, Rcpp::String elem) {
  bool res = std::find(vec.begin(), vec.end(), elem) != vec.end();
  return(res);
}

#endif