#include <Rcpp.h>
#include "helpers.h"
#include "interface.h"

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".mackBoot")]]
Rcpp::NumericVector mack_boot(Rcpp::NumericMatrix triangle, int n_boot, int resids_type, int boot_type, int dist) {
    int n_dev = triangle.rows();
    double** triangle_f = array_create(n_dev, n_dev);

    for (int i = 0; i < n_dev; i++) {
      for (int j = 0; j < n_dev; j++) {
        triangle_f[j][i] = triangle(i, j);
      }
    }

    Rcpp::NumericVector reserve(n_boot);
    mack_boot_f(n_dev, &triangle_f[0][0], resids_type, boot_type, dist, n_boot, reserve.begin());
    free(triangle_f);
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
// [[Rcpp::export(name = ".mackSim")]]
Rcpp::NumericMatrix mack_sim(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::NumericMatrix config, int type) {

  int n_dev = triangle.rows();

  int n_config = config.rows();
  int m_config = config.cols();

  double** triangle_f = array_create(n_dev, n_dev);

  for (int i = 0; i < n_dev; i++) {
    for (int j = 0; j < n_dev; j++) {
      triangle_f[j][i] = triangle(i, j);
    }
  }

  double** config_f = array_create(m_config, n_config);

  for (int i=0; i < n_config; i++) {
    for (int j=0; j < m_config; j++) {
      config_f[j][i] = config(i, j);
    }
  }

  int n_results = n_boot * n_config;
  int m_results = m_config + 1;

  double** results_f = array_create(m_results, n_results);

  mack_sim_f(n_dev, &triangle_f[0][0], n_config, m_config, &config_f[0][0], type, n_boot, &results_f[0][0]);

  Rcpp::NumericMatrix results(n_results, m_results);

    for (int i=0; i < n_results; i++) {
      for (int j=0; j < m_results; j++) {
        results(i, j) = results_f[j][i];
      }
    }


  free(triangle_f);
  free(config_f);
  free(results_f);

  return(results);
};
