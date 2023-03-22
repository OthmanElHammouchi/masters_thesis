#include <Rcpp.h>
#include "helpers.h"

extern "C" {
  
  void glm_boot_f(int n_dev, double* triangle, int n_boot, double* reserve);

  void glm_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  int validate_rng_f(int n_samples);
}

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".glmBoot")]]
Rcpp::NumericVector glm_boot(Rcpp::NumericMatrix triangle, int n_boot) {
    int n_dev = triangle.rows();
    Rcpp::NumericVector reserve(n_boot);
    glm_boot_f(n_dev, triangle.begin(), n_boot, reserve.begin());
  return(reserve);
};

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".glmSim")]]
Rcpp::NumericMatrix glm_sim(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::NumericMatrix config, int type) {
  int n_dev = triangle.rows();
  int n_config = config.rows();
  int m_config = config.cols();
  int n_results = n_boot * n_config;
  int m_results = m_config + 1;
  Rcpp::NumericMatrix results(n_results, m_results);
  glm_sim_f(n_dev, triangle.begin(), n_config, m_config, config.begin(), type, n_boot, results.begin());
  return(results);
};