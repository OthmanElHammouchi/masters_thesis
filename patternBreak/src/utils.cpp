#include <Rcpp.h>
#include "interface.h"

using namespace Rcpp;

// [[Rcpp::export]]
String validate_rng(int n_samples) {
  int success = validate_rng_f(n_samples);
  if (success == 0) {
    return(String("Success"));
  } else {
    return(String("Failure"));
  }
}