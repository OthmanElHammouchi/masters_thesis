#include <math.h>
#include <Rcpp.h>

#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <trng/uniform01_dist.hpp>
#include <RcppThread.h>

using namespace Rcpp;
using namespace trng;

extern "C" {

  void mack_boot_f(int n_dev, double* triangle, int resids_type, int boot_type, int dist, int n_boot, double* reserve);

  void glm_boot_f(int n_dev, double* triangle, int n_boot, double* reserve);

  void mack_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  void glm_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  void* init_rng(int seed) {
    void* rng = (void*) new yarn2((unsigned long) seed);
    return(rng);
  }

  void* lrng_create(void* rng, int n_threads, int i_thread) {
    yarn2* ptr = static_cast<yarn2*>(rng);
    yarn2* _lrng = new yarn2(*ptr);
    (*_lrng).split((unsigned int) n_threads, (unsigned int) i_thread);
    void* lrng = (void*) &_lrng;
    return(lrng);
  }

  double rnorm_par(void* lrng, double mean, double sd) {
    normal_dist<> dist(mean, sd);
    yarn2* ptr = static_cast<yarn2*>(lrng);
    double sample = dist(*ptr);
    return(sample);
}

  double rgamma_par(void* lrng, double shape, double scale) {
    gamma_dist<> dist(shape, scale);
    yarn2* ptr = static_cast<yarn2*>(lrng);
    double sample = dist(*ptr);
    return(sample);
  }

  int rpois_par(void* lrng, double mean) {
    poisson_dist dist(mean);
    yarn2* ptr = static_cast<yarn2*>(lrng);
    int sample = dist(*ptr);
    return(sample);
  }

  double runif_par(void* lrng) {
    uniform01_dist<> dist;
    yarn2* ptr = static_cast<yarn2*>(lrng);
    double sample = dist(*ptr);
    return(sample);
  }

  void check_user_input(void) {
    RcppThread::checkUserInterrupt();
  }

  void* pgbar_create(int total, int freq) {
    void* ptr = (void*) new RcppThread::ProgressBar(total, freq);
    return(ptr);
  }

  void pgbar_incr(void* progress_bar) {
    (*static_cast<RcppThread::ProgressBar*>(progress_bar))++;
  }
}

double** array_create(int n, int m) {
    double* data = (double*)malloc(n * m * sizeof(double));
    double** array = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) {
      array[i] = &(data[m * i]);
    }
    return(array);
}

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = "mackBoot")]]
NumericVector mack_boot(NumericMatrix triangle, int n_boot, int resids_type, int boot_type, int dist) {


    int n_dev = triangle.rows();
    
    double** triangle_f = array_create(n_dev, n_dev);

    for (int i = 0; i < n_dev; i++) {
      for (int j = 0; j < n_dev; j++) {
        triangle_f[j][i] = triangle(i, j);
      }
    }

    NumericVector reserve(n_boot);

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
// [[Rcpp::export(name = "glmBoot")]]
NumericVector glm_boot(NumericMatrix triangle, int n_boot) {

    int n_dev = triangle.rows();
    
    double** triangle_f = array_create(n_dev, n_dev);

    for (int i = 0; i < n_dev; i++) {
      for (int j = 0; j < n_dev; j++) {
        triangle_f[j][i] = triangle(i, j);
      }
    }

    NumericVector reserve(n_boot);

    glm_boot_f(n_dev, &triangle_f[0][0], n_boot, reserve.begin());

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
NumericMatrix mack_sim(NumericMatrix triangle, int n_boot, NumericMatrix config, int type) {

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

  NumericMatrix results(n_results, m_results);

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

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".glmSim")]]
NumericMatrix glm_sim(NumericMatrix triangle, int n_boot, NumericMatrix config, int type) {

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

  glm_sim_f(n_dev, &triangle_f[0][0], n_config, m_config, &config_f[0][0], type, n_boot, &results_f[0][0]);

  NumericMatrix results(n_results, m_results);

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

// [[Rcpp::export]]
NumericVector test() {
  NumericVector result(1e5);
  yarn2 rng(42);
  for (int i; i < 1e5; i++) {
    poisson_dist dist(i);
    result(i) = dist(rng);
  };
  return(result);
};

// [[Rcpp::export]]
NumericVector test2() {
  NumericVector result(1e5);
  yarn2 rng(42);
  for (int i; i < 1e5; i++) {
    normal_dist<> dist(i, i);
    result(i) = dist(rng);
  };
  return(result);
};