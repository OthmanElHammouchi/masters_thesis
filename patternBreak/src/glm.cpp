#include <Rcpp.h>
#include "helpers.h"

extern "C" {
  
  void glm_boot_(int n_dev, double* triangle, int n_boot, double* reserve);

  void glm_sim_(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  int validate_rng_(int n_samples);
}

//' Simulate Poisson GLM reserve.
//'
//' @param triangle Incremental claims triangle
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector glm_boot(Rcpp::NumericMatrix triangle, int n_boot) {
    int n_dev = triangle.rows();
    Rcpp::NumericVector reserve(n_boot);
    glm_boot_(n_dev, triangle.begin(), n_boot, reserve.begin());
  return(reserve);
};

//' Simulate Poisson GLM reserve for different perturbed and excluded points.
//'
//' @param triangle Incremental claims triangle
//' @param sim_type Simulation type: `"single"` (the default), `"calender"`, or `"origin"`
//' @param n_boot Number of bootstrap simulations
//' @param factor Perturbation factor
//' @param boot_type Type of bootstrap: `"parametric"`, `"residuals"`, or `"pairs"`
//' @param proc_dist Distribution of process error: `"normal"` or `"gamma"`
//' @param cond Specified whether the bootstrap should be conditional or unconditional. Default is `TRUE`.
//'
//' @details
//' The simulation configuration inputs `sim_type`, `factor`, `boot_type`, `proc_dist` and `conditional` can be
//' either strings or character vectors. In the latter case, the simulation is computed for all feasible combinations.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix glm_sim(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::NumericMatrix config, int type) {
  int n_dev = triangle.rows();
  int n_config = config.rows();
  int m_config = config.cols();
  int n_results = n_boot * n_config;
  int m_results = m_config + 1;
  Rcpp::NumericMatrix results(n_results, m_results);
  glm_sim_(n_dev, triangle.begin(), n_config, m_config, config.begin(), type, n_boot, results.begin());
  return(results);
};