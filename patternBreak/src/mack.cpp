#include <Rcpp.h>
#include "helpers.h"
#include "interface.h"
#include "constants.h"

//' Simulate Mack CL reserve.
//'
//' @param triangle Cumulative claims triangle
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mackBoot(
  Rcpp::NumericMatrix triangle,
  int n_boot, 
  Rcpp::String boot_type, 
  Rcpp::String process_dist, 
  bool conditional, 
  Rcpp::String resids_type) {
  int n_dev = triangle.rows();
  double ** triangle__ = array_from_rcpp(triangle);
  double** triangle_ = array_transpose(n_dev, n_dev, triangle__);

  double* reserve_ = (double*) malloc(sizeof(double) * n_boot); 
  mack_boot_(n_dev, &triangle_[0][0], key[boot_type], key[process_dist], conditional, key[resids_type], n_boot, &reserve_[0]);

  Rcpp::NumericVector reserve(n_boot);
  for (int i=0; i < n_boot; i++) {
    reserve(i) = reserve_[i];
  }

  free(triangle_);
  free(reserve_);

  return(reserve);
};

//' Simulate Mack CL reserve for different perturbed and excluded points.
//'
//' @param triangle Cumulative claims triangle
//' @param sim_type Simulation type: `"single"` (the default), `"calender"`, or `"origin"`
//' @param n_boot Number of bootstrap simulations
//' @param factors Vector of perturbation factors
//' @param boot_type Type of bootstrap: `"parametric"`, `"non-parametric"`, or `"pairs"`
//' @param process_dist Distribution of process error: `"normal"` or `"gamma"`
//' @param conditional Specified whether the bootstrap should be conditional or unconditional. Default is `TRUE`.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix mackSim(
  Rcpp::NumericMatrix triangle, 
  Rcpp::String sim_type,
  int n_boot, 
  Rcpp::NumericVector factors,
  Rcpp::String boot_type,
  Rcpp::String process_dist,
  bool conditional = true,
  Rcpp::String resids_type = "none",
  bool show_progress = true
  ) {
  int n_dev = triangle.rows();
  double** triangle_ = array_transpose(n_dev, n_dev, array_from_rcpp(triangle));

  int n_res, m_res;
  if (key[sim_type] == SINGLE) {
    int n_outliers = (pow(n_dev, 2) - n_dev) / 2 - 1;
    int n_excl = n_outliers;
    int n_factors = factors.length();
    n_res = n_boot * n_outliers * n_excl * n_factors;
    m_res = 6;
  } else if (key[sim_type] == CALENDAR) {
    int n_outliers = n_dev - 1;
    int n_excl = n_outliers;
    int n_factors = factors.length();
    n_res = n_boot * n_outliers * n_excl * n_factors;
    m_res = 4;
  } else if (key[sim_type] == ORIGIN) {
    int n_outliers = n_dev;
    int n_excl = n_outliers;
    int n_factors = factors.length();
    n_res = n_boot * n_outliers * n_excl * n_factors;
    m_res = 4;
  }

  double** sim__ = array_create(m_res, n_res);
  mack_sim_(n_dev, &triangle_[0][0], key[sim_type], n_boot, factors.length(), factors.begin(), key[boot_type], key[process_dist], conditional, key[resids_type], show_progress, n_res, m_res, &sim__[0][0]);

  double** sim_ = array_transpose(m_res, n_res, sim__);

  Rcpp::NumericMatrix sim(n_res, m_res);

  for (int i=0; i < n_res; i++) {
    for (int j=0; j < m_res; j++) {
      sim(i, j) = sim_[i][j];
    }
  }

  free(triangle_);
  free(sim_);
  free(sim__);

  return(sim);
};
